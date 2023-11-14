
# setup -------------------------------------------------------------------

pacman::p_load(tidyverse,
               lme4, stats)

# Poisson GLMM --------------------------------------------------------------------

# data
set.seed(123)

n_sample <- 100
size <- 10

x1 <- rpois(n_sample, 1)
x2 <- rbinom(n_sample, 1, prob = 0.5)
x3 <- sample(x = c("site1", "site2", "site3", "site4", "site5"),
             prob = c(.1, .2, .1, .3, .3),
             size = n_sample,
             replace = T)
x3_num <- as.numeric(factor(x3))

# parameter values

## b = regression coef
## g = group level errors (i.e., random effects)
b <- c(0.1, 0.2, 0.3)
eps <- rnorm(n_distinct(x3_num), mean = 0, sd = 1)

## simulate data
X <- model.matrix(~ x1 + x2)
y_hat <- boot::inv.logit(X %*% b + eps[x3_num]) # needs to be binary
y <- rbinom(n = nrow(X), size = size , prob = y_hat)  # needs to be positive but how 

data1 <- data.frame(x1 = x1,
                    x2 = x2,
                    x3 = x3,
                    y = y,
                    size = size)

# Generalized linear mixed model 
glmer(cbind(y, size - y) ~ x1 + x2 + (1 | x3), # random intercept
      data = data1,
      family = "binomial") %>% 
  summary()

## jags model setup for glmm--------------------------------------------------------------

## parameters ####
para <- c("b0",
          "b1",
          "b2",
          "tau1",
          "eps")

## model file ####
m3 <- runjags::read.jagsfile("code/example_model_binomial.R")

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


## jags for glmm--------------------------------------------------------------------

d_jags3 <- list(Y = y,
                Size = size,
                X1 = x1,
                X2 = x2,
                G = x3_num,
                Nsample = length(y),
                Ngroup = n_distinct(x3_num))

## run jags ####
post3 <- runjags::run.jags(model = m3$model,
                           monitor = para,
                           data = d_jags3,
                           n.chains = n_chain,
                           inits = inits,
                           method = "parallel",
                           burnin = n_burn,
                           sample = n_sample,
                           adapt = n_ad,
                           thin = n_thin,
                           n.sims = n_chain,
                           module = "glm") #specific to jags doesnt need to be change

(mcmc_summary3 <- MCMCvis::MCMCsummary(post3$mcmc))

