#example_model_binomial

# Model for binomial GLMM ----------------------------------------------------------

model {
  
  ninfo <- 0.01
  
  # priors ------------------------------------------------------------------
  
  ## precision
  tau1 ~ dgamma(0.01, 0.01) 
  
  b0 ~ dnorm(0, ninfo)
  b1 ~ dnorm(0, ninfo)
  b2 ~ dnorm(0, ninfo)
  
  
  # likelihood --------------------------------------------------------------
  
  ## for data-level replicates
  for (i in 1:Nsample) {
    
    Y[i] ~ dbin(p[i], Size)
    #Y[i] ~ dbern(p[i], y[i]) ## need probability term but which is this on run model?
    logit(p[i]) <- logit.p[i]  #link function
    logit.p[i] <- (b0 + eps[G[i]]) + b1 * X1[i] + b2 * X2[i]     #linear predictor and random effect
  }
  
  ## for group-level replicates (random effects)
  ## tau1 is the precision for the random effect
  for (j in 1:Ngroup) {       # define group effect
    eps[j] ~ dnorm(0, tau1)
  }
  
}