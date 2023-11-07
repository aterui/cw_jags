# Model for poisson GLMM ----------------------------------------------------------

model {
  
  ninfo <- 0.01
    
  # priors ------------------------------------------------------------------
  
  ## precision
  tau1 ~ dgamma(0.01, 0.01) 
    
  ## precision to SDs
  sd1 <- sqrt(1 / tau1)
  
  b0 ~ dnorm(0, ninfo)
  b1 ~ dnorm(0, ninfo)
  b2 ~ dnorm(0, ninfo)
  
     
  # likelihood --------------------------------------------------------------
  
  ## for data-level replicates
  for (i in 1:Nsample) {
    
    Y[i] ~ dpois(lambda[i]) #lambda here is the term called for mean of poisson ie distribution
    log(lambda[i]) <- log.lambda[i]  #link function
    log.lambda[i] <- (b0 + eps[G[i]]) + b1 * X1[i] + b2 * X2[i]     #linear predictor and random effect
  }
  
  ## for group-level replicates (random effects)
  ## tau1 is the precision for the random effect
  for (j in 1:Ngroup) {       # define group effect
    eps[j] ~ dnorm(0, tau1)
  }
    
}

