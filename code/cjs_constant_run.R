# real cjs constant model run
# setup -------------------------------------------------------------------

pacman::p_load(tidyverse,
               lme4)

source("code/real_function.R")
get_first <- function(x) min(which(!is.na(x)))

# cjs simulated data ------------------------------------------------------

## see `code/real_function.R` for source codes

Y <- df_0$recap
Nind <- nrow(df_0)
Nocc <- (df_0$occasion)
Fc <- apply(df_0, MARGIN = 1, get_first)

d_jags <- list(Y = df_0$recap,       # recapture 0-1 
               Nind = Nind,          # number of individuals
               Nocc = Nocc, # occasion 1-13
               Fc = Fc)              # first occasion 

para <- c("mean.phi", "mean.p")

# mcmc setup --------------------------------------------------------------

## model file ####
mcjs <- runjags::read.jagsfile("code/cjs_constant_model.R")

## mcmc setup ####
n_ad <- 1000
n_iter <- 2.0E+3
n_thin <- max(3, ceiling(n_iter / 250))
n_burn <- ceiling(max(10, n_iter/2))
n_sample <- ceiling(n_iter / n_thin)
n_chain <- 4

inits <- replicate(n_chain,
                   list(.RNG.name = "base::Mersenne-Twister",
                        .RNG.seed = NA,
                        mean.p = 0.5,
                        mean.phi = 0.9),
                   simplify = FALSE)

for (j in 1:n_chain) inits[[j]]$.RNG.seed <- (j - 1) * 10 + 1


# run ---------------------------------------------------------------------

post <- runjags::run.jags(model = mcjs$model,
                          monitor = para,
                          data = d_jags,
                          n.chains = n_chain,
                          inits = inits,
                          method = "parallel",
                          burnin = n_burn,
                          sample = n_sample,
                          adapt = n_ad,
                          thin = n_thin,
                          n.sims = n_chain,
                          module = "glm") #specific to jags doesnt need to be change

MCMCvis::MCMCsummary(post$mcmc)
