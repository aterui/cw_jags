# Model for GLMM ----------------------------------------------------------

model {
  
  # priors ------------------------------------------------------------------
  
  ## precision
  tau0 ~ dgamma(0.01, 0.01) # tau here is a precision parameter that cannot be negative --> errors will occur if negative
  tau1 ~ dgamma(0.01, 0.01) # tau here is a precision parameter that cannot be negative --> errors will occur if negative
  
  ## precision to SDs
  sd0 <- sqrt(1 / tau0)
  sd1 <- sqrt(1 / tau1)
  
  b0 ~ dnorm(0, 0.01)
  b1 ~ dnorm(0, 0.01)
  b2 ~ dnorm(0, 0.01)
  # a <- 1 / (sd*sd)
  # sd ~ dnorm(0,5)
  
    
  # likelihood --------------------------------------------------------------
  
  ## for data-level replicates
  for (i in 1:Nsample) {
    Y[i] ~ dnorm(mu[i], tau0)
    mu[i] <- b0 + b1 * X1[i] + b2 * X2[i] + eps[G[i]]
    
    # Y[i] ~ dpois(lambda[i]) # lambda here is the term called for mean of poisson
    # log(lambda[i]) <- log.lambda[i]
    # log.lambda[i] <- b0 + b1 * X1[i] + b2 * X2[i] + eps[i]
    # eps[i] ~ dnorm(0, a)
  }
  
  ## for group-level replicates (random effects)
  ## tau1 is the precision for the random effect
  for (j in 1:Ngroup) {
    eps[j] ~ dnorm(0, tau1)
  }
    
}

