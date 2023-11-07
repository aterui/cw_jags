# Model for poisson GLMM ----------------------------------------------------------

model {
  
  # priors ------------------------------------------------------------------
  
  ## precision
  tau0 ~ dgamma(0.01, 0.01) # tau here is a precision parameter that cannot be negative 
  tau1 ~ dgamma(0.01, 0.01) 
 
    
  ## precision to SDs
  sd0 <- sqrt(1 / tau0)   # allows for better transformation / understanding of model output 
  sd1 <- sqrt(1 / tau1)
  
  b0 ~ dpois(100, 1) 
  b1 ~ dpois(100, 1)
  b2 ~ dpois(100, 1)
  
     
  # likelihood --------------------------------------------------------------
  
  ## for data-level replicates
  for (i in 1:Nsample) {
    
    Y[i] ~ dpois(lambda[i], tau0) #lambda here is the term called for mean of poisson ie distribution
    log(lambda[i]) <- exp(log.lambda[i])  #link function
    log.lambda[i] <- b0 + b1 * X1[i] + b2 * X2[i] + eps[G[i]]     #linear predictor and random effect
  }
  
  ## for group-level replicates (random effects)
  ## tau1 is the precision for the random effect
  for (j in 1:Ngroup) {       # define group effect
    eps[j] ~ dnorm(0, tau1)
  }
    
}

