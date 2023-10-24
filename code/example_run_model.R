# GLMM --------------------------------------------------------------------

# data
set.seed(123)
x1 <- rnorm(1000, 0, 1)
x2 <- rnorm(1000, 0, 1)
x3 <- sample(x = c("site1", "site2", "site3", "site4", "site5"),
             prob = c(.1, .2, .1, .3, .3),
             size = 100,
             replace = T)

b <- c(0.1, 0.2, 0.3)

X <- model.matrix(~ x1 + x2)
y_hat <- X %*% b
y <- rnorm(n = nrow(X), mean = y_hat, sd = 1)


data1 <- cbind(x1, x2, x3, y) %>% 
  data.frame()  

data1$x1 <- as.numeric(as.character(data1$x1))
data1$x2 <- as.numeric(as.character(data1$x2))
data1$y <- as.numeric(as.character(data1$y))

# Generalized linear mixed model 
lmer(y ~ x1 + x2 + (1 | x3), # random intercept
     data = data1) %>% 
  summary()

# group = random effect term
# (1 | group) = random intercept
# (0 + x | group) or (-1 + x | group) = random slope
# (x | group) or (1 + x | group ) = random slope and correlated intercept
# (1 | group) or (0 + x | group) = random slope and uncorrelated intercept


## jags model setup for glmm--------------------------------------------------------------

## parameters ####
para <- c("tau",
          "b0",
          "b1",
          "b2",
          "lambda",
          "sd",
          "eps")

## model file ####
m3 <- runjags::read.jagsfile("example_model.R")

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
                X1 = x1,
                X2 = x2,
                Nsample = length(y))

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

mcmc_summary3 <- MCMCvis::MCMCsummary(post3$mcmc)

