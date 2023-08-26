library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(GGally)
library(gridExtra)
library(VGAM)
library(geometry)

# Importing data
d <- read_delim('data/winequality-white.csv', delim = ";")
y <- d$quality
X <- d[, -12]
X <- matrix(unlist(X), ncol = ncol(X), byrow = FALSE)
X <- cbind(rep(1,nrow(X)),scale(X))
X_names <- c("b0", colnames(d)[1:11])
colnames(X) <- X_names

# Plotting histogram of quality scores
y_summary <- as_tibble(y)
hist_plot <- ggplot(y_summary) + geom_bar(aes(x=y), fill = "#4B7599", color = "#4B7599") +
  ggtitle('White wine quality scores histogram')+ theme_classic() +
  scale_x_continuous(breaks = seq(y[which.min(y)],y[which.max(y)],1))
hist_plot

# Saving it
ggsave(filename = paste('plots/White wine quality scores histogram.png'),
       plot = hist_plot, width = 7, height = 7, limitsize = TRUE)

# Plotting a correlation matrix
correlation_plot <- ggcorr(X[,2:ncol(X)], palette = "RdYlGn", name = "rho",
       label = TRUE, label_color = "black",
       label_alpha = TRUE) + ggtitle('White wine. Correlation of explanatory variables')
correlation_plot

# Saving it
ggsave(filename = paste('plots/White wine correlation matrix.png'),
       plot = correlation_plot, width = 10, height = 10, limitsize = TRUE)


n <- length(y)
k <- ncol(X)
npars <- k + (k-1) + 2

logprior <- function(pars) {
  #' The functions calculates log probability density function
  #' of beta (Laplace distribution) and sigma2 (Inverse Gamma) priors
  #' It should be noted that Laplace destribution of beta priors is
  #' represented by hiararchical model that depends on parameters
  #' lambda2 (Gamma distributed), tau (Exponencially distributed)
  #' and sigma2 (distributed as mentioned above)
  b <- as.matrix(pars[1:k])
  tau <- pars[(k+1):(2*k-1)]
  lambda2 <- pars[2*k]
  sigma2 <- pars[npars]
  if (lambda2 <= 0)
    return(log(0))
  if (sigma2 <= 0)
    return(log(0))
  if (any(tau<0))
    return(log(0))
  dens <- (
    sum(dnorm(abs(b[2:k]),0,sqrt(sigma2*tau),log=TRUE))+
      dgamma(1/sigma2,4,60,log=TRUE) +
      dgamma(lambda2,0.025,0.1,log=TRUE)+
      sum(dexp(tau,lambda2/2,log=TRUE)))
  return(dens)
}

loglikelihood <- function(y,X,par) {
  #'    Calculation of likelihood function for response variable y assuming
  #'    linear regression with y distributed Normally with
  #'    mean equal to sum of predictors X multiplied by regression
  #'    coefficients b and standard deviation sigma
  b <- as.matrix(par[1:k])
  sigma2 <- par[npars]
  if (sigma2 <= 0)
    return(log(0))
  dens <- sum(dnorm(y,t(X%*%b),sigma2,log=TRUE))
  return(dens)
}

# Preparation of matris of draws
number_of_draws <- 120000
draws <- matrix(0,nrow=number_of_draws,ncol=npars)
draws[1,1] <- mean(y) # first beta is intercept parameter that is initially consider equalt mean od response variable
draws[1,2:k] <- rep(0,k-1)  # beta coefficients of explanatory variables
draws[1,(k+1):(2*k-1)] <- rep(1,k-1) # tau - parameter in Laplace distribution (see project report)
draws[1,npars-1] <- 1 # lambda2 - parameter in Laplace distribution (see project report)
draws[1,npars] <- 1/rgamma(1,4,60) # sigma2

n_accepted <- 0 # to check how many steps in sampling below were accepted

for (step in 2:number_of_draws) {
  proposed <- rnorm(npars,draws[step-1,],0.01)
  r <- loglikelihood(y, X, proposed) +
    logprior(proposed) -
    loglikelihood(y,X,draws[step-1,]) -
    logprior(draws[step-1,])
  u <- runif(1)
  if (log(u) < r) {
    draws[step,] <- proposed
    n_accepted <- n_accepted+1
  } else {
    draws[step,] <- draws[step-1,]
  }
}


# To load the draws that I used in the report (one out of the three below, not all at the same time)
# load("draws/draws_white_normal 1.RData")
# load("draws/draws_white_normal 2.RData")
# load("draws/draws_white_normal 3.RData")

sample <- draws[50000:number_of_draws,] # burn-in tail of first 50,000 samples

# Preparation of posterior ditribution of betas for plotting
betas <- sample[,2:k]
colnames(betas) <- X_names[2:k]
betas <- as_tibble(betas)
betas_gathered <- gather(betas, key = 'beta_name', value = 'beta_value')

# Plotting poterior distribution of betas and saving the plot
beta_posteriors_plot <- ggplot(betas_gathered, aes(beta_value, beta_name)) +
  geom_violin(aes(), fill = "#4B7599", color = "#4B7599") +
  ggtitle('Posterior distribution of betas for white wine with Normal response model')
beta_posteriors_plot

ggsave(filename = paste('plots/White - Normal. Posterior distribution of betas Draw #3.png'),
       plot = beta_posteriors_plot, height = 10, width = 8)

# Sample from posterior distribution
# First creating possible values of wine quality
# As in the using linear model to implement the inference
# We consider that 'y' = quality of wine is continuous value
# in the range from 0 to 10.
# Therefore we create a grid of possible values of y variable:
ygrid <- seq(from=0, to=10, by=0.1)
ydensity <- rep(0, length(ygrid))
ppd_tables <- list() # to collect posterior predictive disrtribution for selected wines in table format
ppd_plots <- list() # to collect plots posterior predictive disrtibution for selected wines

posterior <- sample*0 # same size matrix as sample of posterior but filled with 0s
wine_check <- c(1178, 2374, 4619)

# Plotting posterior predictive distribution of 3 wines above
# Loops below will take a while to finish
for(wine in seq_along(wine_check)) {
  for(i in seq_along(ygrid)) {
        for(j in seq_len(nrow(posterior))) {
        # integrating over all posterior
        # distribution of beta coefficients
          posterior[j,] <- exp(loglikelihood(ygrid[i], X[wine_check[wine],], sample[j,]))
        }
    ydensity[i] <- sum(posterior)
  }
  ydensity <- ydensity/sum(ydensity) # normalizing posterior predictive distribution
  ppd_tables[[wine]] <- tibble(wine_quality = ygrid, probability = ydensity)
  ppd_plots[[wine]] <- ggplot(ppd_tables[[wine]], aes(x=wine_quality, y=probability)) +
    geom_bar(stat = 'identity', fill = '#4B6A99') +
    geom_vline(xintercept = y[wine_check[wine]], size = 2, color = '#AB200E') +
    scale_x_continuous(breaks = seq(0,10,1)) +
    ggtitle(paste("Posterior predictive distribution (normal response - original variables) for white wine number",
                  wine_check[wine]))
}

# Showing the plots generated above and saving it
ppd_plots[[1]]
ppd_plots[[2]]
ppd_plots[[3]]
for(i in seq_along(ppd_plots)) ggsave(filename = paste('plots/Normal response and original variables',
                                                       'PPD of white wine #',wine_check[i],'.png'),
                                      plot = ppd_plots[[i]])

# Transformation of variables
# TRANSFORMED VARIABLE 1 - Proportion of free free sulphur dioxide in total sulfur dioxide
X_plot <- as_tibble(X)
sulfur_dioxide_plot <- qplot(x = `free sulfur dioxide`,y = `total sulfur dioxide`, data=X_plot) +
  geom_rug(col=rgb(.5,0,0,alpha=.2))
sulfur_dioxide_plot

ggsave('plots/White wine. free and total sulphur dioxide.png', plot = sulfur_dioxide_plot)

# TRANSFORMED VARIABLE 2 - Volatile / Fixed acidity
vol_fixed_acidity_plot <- qplot(x = `volatile acidity`,y = `fixed acidity`, data=X_plot) +
  geom_rug(col=rgb(.5,0,0,alpha=.2))
vol_fixed_acidity_plot

ggsave('plots/White wine. Volatile and fixed acidity.png', plot = vol_fixed_acidity_plot)

# TRANSFORMED VARIABLE 3 - Log transformation of chlorides
chlorides_plot <- ggplot(X_plot, aes(x=`chlorides`, y=y)) +
  geom_point() + geom_smooth(method=lm)
chlorides_plot

ggsave('plots/White wine. Chlorides and quality scatterplot.png', plot = chlorides_plot)

log_chlorides_plot <- ggplot(X_plot, aes(x=log(`chlorides`), y=y)) +
  geom_point() + geom_smooth(method=lm)
log_chlorides_plot

ggsave('plots/White wine. Log Chlorides and quality scatterplot.png', plot = log_chlorides_plot)

# Modifying our variables dataframe:
X_upd <- d[, -12]

X_upd <- X_upd %>% mutate(`free to total SD` = `free sulfur dioxide` / `total sulfur dioxide`,
                          `volatile to fixed acidity` = `volatile acidity` / `fixed acidity`,
                          `log chlorides` = log(`chlorides`)) %>%
  select(-`free sulfur dioxide`, -`total sulfur dioxide`, -`volatile acidity`, -`fixed acidity`, -`chlorides`)

X_upd_names <- c("b0", colnames(X_upd))
X_upd <- cbind(rep(1,nrow(X_upd)),scale(X_upd))
colnames(X_upd) <- X_upd_names
as_tibble(X_upd)
correlation_plot_upd <- ggcorr(X_upd[,2:ncol(X_upd)], palette = "RdYlGn", name = "rho",
       label = TRUE, label_color = "black",
       label_alpha = TRUE) + ggtitle('White wine. Correlation of explanatory variables (transformed)')
correlation_plot_upd

ggsave(filename = paste('plots/White wine correlation matrix (transformed variables).png'),
       plot = correlation_plot_upd, width = 10, height = 10, limitsize = TRUE)

# Starting from here, the code is the same as above code for original variables
# Doing the same that we performed with original variables
n <- length(y)
k <- ncol(X_upd)
npars <- k + (k-1) + 2

number_of_draws <- 120000
draws <- matrix(0,nrow=number_of_draws,ncol=npars)
draws[1,1] <- mean(y)
draws[1,2:k] <- rep(0,k-1)
draws[1,(k+1):(2*k-1)] <- rep(1,k-1)
draws[1,npars-1] <- 1
draws[1,npars] <- 1/rgamma(1,4,60)

n_accepted <- 0

for (step in 2:number_of_draws) {
  proposed <- rnorm(npars,draws[step-1,],0.01)
  r <- loglikelihood(y, X_upd, proposed) +
    logprior(proposed) -
    loglikelihood(y,X_upd,draws[step-1,]) -
    logprior(draws[step-1,])
  u <- runif(1)
  if (log(u) < r) {
    draws[step,] <- proposed
    n_accepted <- n_accepted+1
  } else {
    draws[step,] <- draws[step-1,]
  }
}

# To load draws used in the report
# load("draws/draws_white_normal_transformed_vars.RData")
sample_upd <- draws[50000:number_of_draws,]
betas_upd <- sample_upd[,2:k]
colnames(betas_upd) <- X_upd_names[2:k]

betas_upd <- as_tibble(betas_upd)
betas_gathered_upd <- gather(betas_upd, key = 'beta_name', value = 'beta_value')
beta_posteriors_plot_upd <- ggplot(betas_gathered_upd, aes(beta_value, beta_name)) +
  geom_violin(aes(), fill = "#4B7599", color = "#4B7599")
beta_posteriors_plot_upd

ggsave(filename = paste('plots/White - Normal. Posterior distribution of betas (transformed variables).png'),
       plot = beta_posteriors_plot_upd)


# Sample from posterior distribution
# First creating possible values of wine quality
# As in the using linear model to implement the inference
# We consider that 'y' = quality of wine is continuous value
# in the range from 0 to 10.
# Therefore we create a grid of possible values of y variable:
ygrid <- seq(from=0, to=10, by=0.1)
ydensity <- rep(0, length(ygrid))
ppd_tables_upd <- list() # to collect posterior predictive disrtribution for selected wines in table format
ppd_plots_upd <- list() # to collect plots posterior predictive disrtibution for selected wines

posterior <- sample_upd*0 # same size matrix as sample of posterior but filled with 0s
wine_check <- c(1178, 2374, 4619)

# Loops below will take a while
for(wine in seq_along(wine_check)) {
  for(i in seq_along(ygrid)) {
        for(j in seq_len(nrow(posterior))) {
        # integrating over all posterior
        # distribution of beta coefficients
          posterior[j,] <- exp(loglikelihood(ygrid[i], X_upd[wine_check[wine],], sample_upd[j,]))
        }
    ydensity[i] <- sum(posterior)
  }
  ydensity <- ydensity/sum(ydensity) # normalizing posterior predictive distribution
  ppd_tables_upd[[wine]] <- tibble(wine_quality = ygrid, probability = ydensity)
  ppd_plots_upd[[wine]] <- ggplot(ppd_tables_upd[[wine]], aes(x=wine_quality, y=probability)) +
    geom_bar(stat = 'identity', fill = '#4B6A99') +
    geom_vline(xintercept = y[wine_check[wine]], size = 2, color = '#AB200E') +
    scale_x_continuous(breaks = seq(0,10,1)) +
    ggtitle(paste("Posterior predictive distribution (normal response - transformed variables) for white wine number",
                  wine_check[wine]))
}

# Plotting posterior predictive distributions
ppd_plots_upd[[1]]
ppd_plots_upd[[2]]
ppd_plots_upd[[3]]

# Saving the plots
for(i in seq_along(ppd_plots_upd)) ggsave(filename = paste('plots/Normal response and transformed variables',
                                                           'PPD of white wine #',wine_check[i],'.png'),
                                          plot = ppd_plots_upd[[i]])