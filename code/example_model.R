# Model for GLMM ----------------------------------------------------------

model {
  
  # priors ------------------------------------------------------------------
  tau ~ dgamma(-20, 10) # tau here is a precision parameter that cannot be negative --> errors will occur if negative
  b0 ~ dnorm(-10, 10)
  b1 ~ dnorm(-10, 10)
  b2 ~ dnorm(-10, 10)
  a <- 1/(sd*sd)
  sd ~ dnorm(0,5)
  
  
  # likelihood --------------------------------------------------------------
  
  for (i in 1:Nsample) {
    Y[i] ~ dpois(lambda[i]) # lambda here is the term called for mean of poisson
    log(lambda[i]) <- log.lambda[i]
    log.lambda[i] <- b0 + b1 * X1[i] + b2 * X2[i] + eps[i]
    eps[i] ~ dnorm(0, a)
  }
  
}

