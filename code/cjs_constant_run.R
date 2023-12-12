
# setup -------------------------------------------------------------------

pacman::p_load(tidyverse,
               lme4)

source("code/function.R")
get_first <- function(x) min(which(!is.na(x)))

# cjs simulated data ------------------------------------------------------

## run `obj <- cjs_data()` to produce data
## this function returns obj$Z and obj$Y
## `Z` is a matrix of true survival state, 1 = alive, 0 = dead
## `Y` is a matrix of capture state, 1 = captured, 0 = not captured
## for both matrices, rows are individuals and columns are occasions
## see `code/function.R` for source codes

## arguments in cjs_data()
## p: detection probability
## s: survival probability
## n0: number of new tagging per occasion
## m_occ: number of occasion

obj <- cjs_data

Y <- obj$Y
Nind <- nrow(Y)
Nocc <- ncol(Y)
Fc <- apply(Y, MARGIN = 1, get_first)

d_jags <- list(Y = Y,
               Nind = Nind,
               Nocc = Nocc,
               Fc = Fc)

para <- c("mean.phi", "mean.p")

# mcmc setup --------------------------------------------------------------

## model file ####
mcjs <- runjags::read.jagsfile("code/example_model_cjs.R")

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
