
# setup -------------------------------------------------------------------

rm(list = ls())
pacman::p_load(tidyverse)

# data --------------------------------------------------------------------

set.seed(123)
x1 <- rnorm(100, 0, 1)
x2 <- rnorm(100, 0, 1)
           
b <- c(0.1, 0.2, 0.3)

X <- model.matrix(~ x1 + x2)
y_hat <- X %*% b
y <- rnorm(n = nrow(X), mean = y_hat, sd = 1)

lm(y ~ x1 + x2) %>% 
  summary()

# jags setup --------------------------------------------------------------

## parameters ####
para <- c("tau",
          "b0",
          "b1",
          "b2")

## model file ####
m <- runjags::read.jagsfile("code/model_lm.R")

## mcmc setup ####
n_ad <- 1000
n_iter <- 1.0E+4
n_thin <- max(3, ceiling(n_iter / 250))
n_burn <- ceiling(max(10, n_iter/2))
n_sample <- ceiling(n_iter / n_thin)
n_chain <- 4

inits <- replicate(n_chain,
                   list(.RNG.name = "base::Mersenne-Twister",
                        .RNG.seed = NA),
                   simplify = FALSE)

for (j in 1:n_chain) inits[[j]]$.RNG.seed <- (j - 1) * 10 + 1


# jags --------------------------------------------------------------------

d_jags <- list(Y = y,
               X1 = x1,
               X2 = x2,
               Nsample = length(y))

## run jags ####
post <- runjags::run.jags(model = m$model,
                          monitor = para,
                          data = d_jags,
                          n.chains = n_chain,
                          inits = inits,
                          method = "parallel",
                          burnin = n_burn,
                          sample = n_sample,
                          adapt = n_ad,
                          thin = n_thin,
                          n.sims = 3,
                          module = "glm")

mcmc_summary <- MCMCvis::MCMCsummary(post$mcmc)
