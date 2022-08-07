library(tidyverse)
library(ggplot2)
library(GGally)
library(gridExtra)
library(VGAM)
library(geometry)
library(comprehenr) # list comprehensions as in Python

# Importing data
d <- read_delim('data/winequality-red.csv', delim = ";")
y <- d$quality
cat <- sort(unique(y))
# transforming y variable to matrix 1599 (number of observations) X 6 (number of quality scores)
# each row contains one 1 (TRUE) and five 0 (FALSE), TRUE represents falling in category
y_transformed <- matrix(nrow = length(y), ncol = length(cat))
for(i in seq_along(y)) {
  y_transformed[i,] <- to_vec(for(j in seq_along(cat)) y[i]==cat[j])
}
colnames(y_transformed) <- paste("Quality score", cat)


X <- d[, -12]
X <- matrix(unlist(X), ncol = ncol(X), byrow = FALSE)
X <- scale(X)
X_names <- colnames(d)[1:11]
colnames(X) <- X_names

n <- length(y)
k <- ncol(X)
npars <- k + 1 # betas + sigma2
categories <- length(unique(y))

draw_paratemeters <- function(draws_list, step) {
  #' For the purpose of making bayesian inference
  #' with multinomial response model, samples for each
  #' quality score are kept in a seperate element of the list (see draws below)
  #' Each list consists of p matrices, each having n rows and k columns,
  #' Where p is number of different outcomes (quality scores), n - number of draws
  #' and k is the number of parameters (betas and sigma2).
  #' This function takes the row of the same index from each of p matrices
  #' The index of the row correspond to previous step in sampling process. The output is then used
  #' in the function "proposed_pars" documented below
  pars <- matrix(ncol=npars, nrow = categories)
      for(i in seq_along(draws_list)) {
    pars[i,] <- draws_list[[i]][step-1, 1:ncol(pars)]
  }
  return(pars)
}

proposed_pars <- function(previous_step_pars) {
  #' The function takes the output of the function "draw_sample" (i.e. sample from previous step)
  #' and draws new samples distributed with Normal distribution with
  #' mean equal to value of previous step sample and fixed standard distribution
  #' It is needed for the Metropolis Hasting sampling process in order to get
  #' Posterior distribution of model parameters
  new_pars <- matrix(ncol = npars, nrow = categories)
  for(i in 1:nrow(previous_step_pars)) {
    new_pars[i,1:npars] <- rnorm(n = npars, mean = previous_step_pars[i,1:npars],sd = 0.01)
  }
  return(new_pars)
}

logprior <- function(pars) {
  #' The function is for calculation of prior distribution of beta coefficients
  #' Under assumption that betas are distributed accordance to Normal distribution with
  #' Mean 0 and standard deviation that is distributed in accordance with
  #' Inverse Gamma distribution.
  betas <- pars[,1:k]
  sigma2 <- pars[,k+1]
  dens <- sum(dnorm(betas, 0, sqrt(sigma2), log = TRUE)) +
              sum(dgamma(1/sigma2,4,60,log=TRUE))
  return(dens)
}

loglikelihood <- function(y, X, pars) {
  #' The function takes a matrix of drawn beta parameters as input
  #' and first calculates probability of the theta of falling into each of the
  #' wine_scores for each observed X as theta_i = exp(X*beta_i)/sum(exp(X*beta_i).
  #' Resulting matrix of thetas with n rows and p colums, where n is the number observations and p is the
  #' number of outcomes (wine scores), serves as the probability of each of n outcomes
  #' in multinomial distribution. It should be noted, that y variable needs to be in form
  #' of n x p matrix where no is the number of opservations and p is the number of outcomes
  #' and each row contains p-1 zeroes and a single 1 representing the outcome of the observation.
  all_thetas <- exp(X %*% t(pars[,1:(npars-1)]))  # softmax fucnction (not normalized)
  diagonal <- diag(1/rowSums(all_thetas)) # for normalization fo the thethas
  all_thetas <- diagonal %*% all_thetas # normalizing thethas so that sum of probabilities in each row was equal 1

  ######### OLD CODE BEGINNING
# The code was previously used to calculate probabilities for each observation (X)
# Replaced by 3 lines above that gives the same result but gets to this result by more straightforward
# operations of linear algebra allowing to calculate it much faster than using loops below.
#  thetas <- rep(0, nrow(pars))
#  all_thetas <- matrix(nrow = n, ncol = length(cat))
#  for (i in 1:n) {
#    for(j in 1:categories) {
#        thetas[j] <- exp(dot(X[i,], pars[j,1:(npars-1)]))
#      }
#    thetas <- thetas/sum(thetas)
#    all_thetas[i,] <- thetas
#  }
  ######### OLD CODE END

  dens <- dmultinom(x=y, prob=all_thetas, log = TRUE)
  return(dens)
}


# I took took 3 samples 10000 draws each and burned in first 5000 draws
# The resulting posterior distribution of betas of each sample is
# wide and differs between each other, however some of the distribution is similar between
# samples and other draws that I took. For example, mode estimate of alcohold is close to -5 for
# Quality Score 5 and close to 1 for quality score 6. It also important to note that all samples make similar
# posterior bredictive distribution, which is relatively accurate.
number_of_draws <- 10000
draws <- list()
for (i in 1:categories) {
  draws[[i]] <- matrix(0,nrow=number_of_draws,ncol=npars)
  draws[[i]][1,npars] <- 1/rgamma(1,4,60)
}

n_accepted <- 0
for (step in 2:number_of_draws) {
  print(step) # just to keep track if the sampling process
  old_pars <- draw_paratemeters(draws, step)
  proposed <- proposed_pars(old_pars)
  r <- loglikelihood(y_transformed, X, proposed) +
    logprior(proposed) -
    loglikelihood(y_transformed,X, old_pars) -
    logprior(old_pars)
  u <- runif(1)
  if (log(u) < r) {
    n_accepted <- n_accepted+1
    for(i in seq_along(draws)) {
      draws[[i]][step, ] <- proposed[i, ]
    }
  } else {
    for(i in seq_along(draws)) {
      draws[[i]][step, ] <- old_pars[i, ]
    }
  }
}

# To load the draws that I used in the report (one out of the three below, not all at the same time)
# load("draws/draws_red_multinomial 1.RData")
# load("draws/draws_red_multinomial 2.RData")
# load("draws/draws_red_multinomial 3.RData")

sample <- list()
for(i in seq_along(draws)) {
  sample[[i]] <- draws[[i]][5001:number_of_draws, ]
}

# Below code is for extracting and plotting posterior distribution of betas.
betas_wide <- list()
betas_gathered <- list()
betas_posterior_plots <- list()
for(i in seq_along(sample)) {
  betas_wide[[i]] <- sample[[i]][,1:k]
  colnames(betas_wide[[i]]) <- X_names
  betas_wide[[i]] <- as_tibble(betas_wide[[i]])
  betas_gathered[[i]] <- gather(betas_wide[[i]], key = 'Beta_name', value = 'Beta_value')
  betas_posterior_plots[[i]] <- ggplot(betas_gathered[[i]], aes(Beta_value, Beta_name)) +
    geom_violin(aes(), fill = "#4B7599", color = "#4B7599") +
    ggtitle(paste("Posterior distribution of beta coefficients for Quality Score",cat[i]))
}
# Plotting posterior distribution of betas
betas_posterior_plots

# Saving the plot
for(i in seq_along(betas_posterior_plots)) ggsave(filename = paste('plots/Red - Multinomial. Posterior distribution,',
                                                                   'of betas for Quality Score',cat[i],' Draw #1.png'),
                                                  plot = betas_posterior_plots[[i]],
                                                  height = 10, width = 8)

# Calculating posterior predictive distribution
# moditfying the function to calculate theta value of a single X values and
# so that to calculate it over posterior distributions of betas below
calc_thetas <- function(X, pars) {
  #' The function calculates thetha value for one particular wine
  #' integrated over the whole posterior distribution of betas
  #' The output is a vector of length p, where p is the number of outcomes (wine scores)
  #' The vector is normalized so that all of its elements sum up to 1 representing all possible
  #' wine score probabilities.
  thetas <- rep(0, nrow(pars))
  for(j in 1:categories) {
        thetas[j] <- exp(dot(X, pars[j,1:(npars-1)]))
      }
  thetas <- thetas/sum(thetas)
  return(thetas)
}

# Here we are estimating probability of getting each quality score to the wine from dataset
# we selected wines #123, 900 and 1001 to demostrate
# We calculating probabilities seperately for each of the abovmentioned bines over the whole distribution of
# betas from the sample calculated by Metropolis method.
# There is no point to additionally draw multinomial desitribution for each of the wine scope as
# thetas already reflect propability of classification of the score and the result will b the same.
wine_check <- c(123, 900, 1001)
ppd_plots <- list()
ppd_thetas <- list()
for(wine in seq_along(wine_check)) {
  posterior_thetas <- rep(0, length(cat))
  for(i in 1:nrow(sample[[1]])) {
    params <- draw_paratemeters(sample, i+1) # becasue draw_params is taking parameters of prior step
    posterior_thetas <- posterior_thetas + calc_thetas(X[wine_check[wine],], params)
  }
  posterior_thetas <- as_tibble(posterior_thetas/nrow(sample[[1]])) # divide n
  posterior_thetas <- posterior_thetas %>% rename(Probability = value) %>%
    mutate(Quality_score = colnames(y_transformed)) %>% select(Quality_score, Probability)
  ppd_thetas[[wine]] <- posterior_thetas
  ppd_plots[[wine]] <- ggplot(data = posterior_thetas) +
    geom_bar(aes(x=Quality_score, y = Probability), stat = 'identity', fill = "#4B7599", color = "#4B7599") +
    geom_vline(xintercept = paste('Quality score',y[wine_check[wine]]), size = 2, color = '#AB200E') +
    ggtitle(paste("Posterior predictive distribution (multinomial response) for red wine number",
                  wine_check[wine]))
}

# Plotting the results
ppd_plots[[1]]
ppd_plots[[2]]
ppd_plots[[3]]

# Saving plots
for(i in seq_along(ppd_plots)) {
  ggsave(filename = paste('plots/PPD of red wine #',wine_check[i],'Multinomial response.png'),
         plot = ppd_plots[[i]])
}
# Probabilities in table format
ppd_thetas

