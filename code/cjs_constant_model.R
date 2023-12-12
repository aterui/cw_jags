model {
  
  # prior -------------------------------------------------------------------
  
  for (i in 1:Nind){
    for (t in Fc[i]:(Nocc - 1)) {
      phi[i,t] <- mean.phi
      p[i,t] <- mean.p
    } #t
  } #i

  mean.phi ~ dunif(0, 1) # Prior for mean survival
  mean.p ~ dunif(0, 1) # Prior for mean recapture
  

  # likelihood --------------------------------------------------------------
  
  for (i in 1:Nind){
    # Define latent state at first capture
    z[i, Fc[i]] <- 1
    
    for (t in (Fc[i] + 1):Nocc){
      
      # State process
      mu1[i,t] <- phi[i, t - 1] * z[i, t - 1]
      z[i, t] ~ dbern(mu1[i,t])
      
      # Observation process
      mu2[i, t] <- p[i, t - 1] * z[i, t]
      Y[i, t] ~ dbern(mu2[i, t])
    } #t
  } #i
}
