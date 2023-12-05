# Chapter 7 CJS Models from Bayesian Population Analysis
# Author: Marc Kery


# Define parameter values
n.occasions <- 6 # Number of capture occasions 
marked <- rep(50, n.occasions - 1 ) # Annual number of newly marked individuals 
phi <- rep(0.65, n.occasions - 1 ) 
p <- rep(0.4, n.occasions - 1) 


# Define matrices with survival and recapture probabilities 
PHI <- matrix(phi, ncol = n.occasions - 1, nrow = sum(marked)) 
P <- matrix(p, ncol = n.occasions - 1, nrow = sum(marked)) 

# Define function to simulate a capture-history (CH) matrix 
simul.cjs <- function(PHI, P, marked){ 
  n.occasions <- dim(PHI)[2] + 1
  
  CH <- matrix(0, ncol = n.occasions, nrow = sum(marked))
  # Define a vector with the occasion of marking
  mark.occ <- rep(1:length(marked), marked[1:length(marked)])
  # Fill the CH matrix
  for (i in 1:sum(marked)){
    CH[i, mark.occ[i]] <- 1 # Write an 1 at the release occasion
    if (mark.occ[i]==n.occasions) next
    for (t in (mark.occ[i]+1):n.occasions){
      # Bernoulli trial: does individual survive occasion?
      sur <- rbinom(1, 1, PHI[i,t-1])
      if (sur==0) break # If dead, move to next individual
      # Bernoulli trial: is individual recaptured?
      rp <- rbinom(1, 1, P[i,t-1])
      if (rp==1) CH[i,t] <- 1
    } #t
  } #i
  return(CH)
  
}

# Execute function
CH <- simul.cjs(PHI, P, marked)

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Specify model in BUGS language
sink("cjs-c-c.bug")
cat("
  model {
# Priors and constraints
  for (i in 1:nind){
  for (t in f[i]:(n.occasions−1)){
  phi[i,t] <- mean.phi
  p[i,t] <- mean.p
  } #t
  } #i
  
  mean.phi ~ dunif(0, 1) # Prior for mean survival
  mean.p ~ dunif(0, 1) # Prior for mean recapture

# Likelihood
  for (i in 1:nind){
# Define latent state at first capture
  z[i,f[i]] <- 1
  for (t in (f[i]+1):n.occasions){

# State process
  z[i,t] ~ dbern(mu1[i,t])
  mu1[i,t] <- phi[i,t-1] * z[i,t-1]

# Observation process
  y[i,t] ~ dbern(mu2[i,t])
  mu2[i,t] <- p[i,t-1] * z[i,t]
  } #t
  } #i
  }
",fill = TRUE)
sink()
# Bundle data
bugs.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions =
                    dim(CH)[2])


# Function to create a matrix of initial values for latent state z
ch.init <- function(ch, f){
  for (i in 1:dim(ch)[1]){ch[i,1:f[i]] <- NA}
  return(ch)
}
# Initial values
inits <- function(){list(z = ch.init(CH, f), mean.phi = runif(1, 0, 1),
                         mean.p = runif(1, 0, 1))}
# Parameters monitored
parameters <- c("mean.phi", "mean.p")
# MCMC settings
ni <- 10000
nt <- 6
nb <- 5000
nc <- 3

# Call WinBUGS from R (BRT 1 min)
cjs.c.c <- bugs(bugs.data, inits, parameters, "cjs-c-c.bug", n.chains =
                  nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE,
                bugs.directory = bugs.dir, working.directory = getwd())

# Summarize posteriors
print(cjs.c.c, digits = 3)


# Inclusion of Latent State Variable --------------------------------------

# Function to create a matrix with information about known latent state z
known.state.cjs <- function(ch){
  state <- ch
  for (i in 1:dim(ch)[1]){
    n1 <- min(which(ch[i,]==1))
    n2 <- max(which(ch[i,]==1))
    state[i,n1:n2] <- 1
    state[i,n1] <- NA
  }
  state[state==0] <- NA
  return(state)
}

# Bundle data
bugs.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2],
                  z = known.state.cjs(CH))
# Function to create a matrix of initial values for latent state z
cjs.init.z <- function(ch,f){
  for (i in 1:dim(ch)[1]){
    if (sum(ch[i,])==1) next
    n2 <- max(which(ch[i,]==1))
    ch[i,f[i]:n2] <- NA
  }
  for (i in 1:dim(ch)[1]){
    ch[i,1:f[i]] <- NA
  }
  return(ch)
}

# Initial values
inits <- function(){list(z = cjs.init.z(CH, f), mean.phi = runif(1, 0, 1),
                         mean.p = runif(1, 0, 1))}
# Parameters monitored
parameters <- c("mean.phi", "mean.p")
# MCMC settings
ni <- 10000
nt <- 6
nb <- 5000
nc <- 3

# Call WinBUGS from R (BRT <1 min)
cjs.c.c <- bugs(bugs.data, inits, parameters, "cjs-c-c.bug", n.chains =
                  nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE,
                bugs.directory = bugs.dir, working.directory = getwd())

# Summarize posteriors
print(cjs.c.c, digits = 3)


# Models with Time Variation ----------------------------------------------
# assume survival and recapture vary independently over time

# Priors and constraints
for (i in 1:nind){
  for (t in f[i]:(n.occasions-1)){
    phi[i,t] <- alpha[t]
    p[i,t] <- beta[t]
  } #t
} #i

for (t in 1:n.occasions-1){
  alpha[t] ~ dunif(0, 1) # Priors for time-spec. survival
  beta[t] ~ dunif(0, 1) # Priors for time-spec. recapture
}

# Define parameter values
n.occasions <- 20 # Number of capture occasions
marked <- rep(30, n.occasions-1) # Annual number of newly marked
individuals
mean.phi <- 0.65
var.phi <- 1 # Temporal variance of survival
p <- rep(0.4, n.occasions-1)

# Determine annual survival probabilities
logit.phi <- rnorm(n.occasions-1, qlogis(mean.phi), var.phi^0.5)
phi <- plogis(logit.phi)

# Define matrices with survival and recapture probabilities
PHI <- matrix(phi, ncol = n.occasions-1, nrow = sum(marked), byrow = TRUE)
P <- matrix(p, ncol = n.occasions-1, nrow = sum(marked))

# Simulate capture-histories
CH <- simul.cjs(PHI, P, marked)

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Specify model in BUGS language
sink("cjs-temp-raneff.bug")
cat("
model {

# Priors and constraints
for (i in 1:nind){
for (t in f[i]:(n.occasions−1)){
logit(phi[i,t]) <- mu + epsilon[t]
p[i,t] <- mean.p
} #t
} #i

for (t in 1:(n.occasions−1)){
epsilon[t] ~ dnorm(0, tau)
}

#mu ~ dnorm(0, 0.001) # Prior for logit of mean survival
#mean.phi <- 1 / (1+exp(−mu)) # Logit transformation
mean.phi ~ dunif(0, 1) # Prior for mean survival
mu <- log(mean.phi / (1−mean.phi)) # Logit transformation
sigma ~ dunif(0, 10) # Prior for standard deviation
tau <- pow(sigma, −2)
sigma2 <- pow(sigma, 2) # Temporal variance
mean.p ~ dunif(0, 1) # Prior for mean recapture

# Likelihood
for (i in 1:nind){

# Define latent state at first capture
z[i,f[i]] <- 1
for (t in (f[i]+1):n.occasions){
# State process
z[i,t] ~ dbern(mu1[i,t])
mu1[i,t] <- phi[i,t−1] * z[i,t−1]
# Observation process
y[i,t] ~ dbern(mu2[i,t])
mu2[i,t] <- p[i,t−1] * z[i,t]
} #t
} #i
}
",fill = TRUE)
sink()

# Bundle data
bugs.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2],
                  z = known.state.cjs(CH))

# Initial values
inits <- function(){list(z = cjs.init.z(CH, f), mean.phi = runif(1, 0, 1),
                         sigma = runif(1, 0, 10), mean.p = runif(1, 0, 1))}

# Parameters monitored
parameters <- c("mean.phi", "mean.p", "sigma2")

# MCMC settings
ni <- 10000
nt <- 6
nb <- 5000
nc <- 3

# Call WinBUGS from R (BRT 17 min)
cjs.ran <- bugs(bugs.data, inits, parameters, "cjs-temp-raneff.bug",
                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE,
                bugs.directory = bugs.dir, working.directory = getwd())

# Summarize posteriors
print(cjs.ran, digits = 3)

# Produce histogram
hist(cjs.ran$sims.list$sigma2, col = "gray", nclass = 35, las = 1,
     xlab = expression(sigma^2), main = "")
abline(v = var.phi, col = "red", lwd = 2)


# Temporal Covariates ------------------------------------------------------

# Define parameter values
n.occasions <- 20 # Number of capture occasions
marked <- rep(15, n.occasions−1) # Annual number of newly marked
individuals
mean.phi <- 0.65
p <- rep(0.4, n.occasions-1)
beta <- -0.3 # Slope of survival-winter
relationship
r.var <- 0.2 # Residual temporal variance

# Draw annual survival probabilities
winter <- rnorm(n.occasions-1, 0, 1^0.5)
logit.phi <- qlogis(mean.phi) + beta*winter + rnorm(n.occasions-1, 0,
                                                    r.var^0.5)
phi <- plogis(logit.phi)

# Define matrices with survival and recapture probabilities
PHI <- matrix(phi, ncol = n.occasions-1, nrow = sum(marked),
              byrow = TRUE)
P <- matrix(p, ncol = n.occasions-1, nrow = sum(marked))

# Simulate capture-histories
CH <- simul.cjs(PHI, P, marked)
# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)
# Specify model in BUGS language
sink("cjs-cov-raneff.bug")
cat("
model {

# Priors and constraints
for (i in 1:nind){
for (t in f[i]:(n.occasions−1)){
logit(phi[i,t]) <- mu + beta*x[t] + epsilon[t]
p[i,t] <- mean.p
} #t
} #i

for (t in 1:(n.occasions−1)){
epsilon[t] ~ dnorm(0, tau)
phi.est[t] <- 1 / (1+exp(−mu-beta*x[t]−epsilon[t])) # Yearly
survival
}

mu ~ dnorm(0, 0.001) # Prior for logit of mean survival
mean.phi <- 1 / (1+exp(−mu)) # Logit transformation
beta ~ dnorm(0, 0.001)I(−10, 10) # Prior for slope parameter
sigma ~ dunif(0, 10) # Prior on standard deviation
tau <- pow(sigma, −2)
sigma2 <- pow(sigma, 2) # Residual temporal variance
mean.p ~ dunif(0, 1) # Prior for mean recapture

# Likelihood
for (i in 1:nind){

# Define latent state at first capture
z[i,f[i]] <- 1
for (t in (f[i]+1):n.occasions){

# State process
z[i,t] ~ dbern(mu1[i,t])
mu1[i,t] <- phi[i,t−1] * z[i,t−1]

# Observation process
y[i,t] ~ dbern(mu2[i,t])
mu2[i,t] <- p[i,t−1] * z[i,t]
} #t
} #i
}
",fill = TRUE)
sink()
# Bundle data
bugs.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2],
                  z = known.state.cjs(CH), x = winter)


# Initial values
inits <- function(){list(z = cjs.init.z(CH, f), mu = rnorm(1), sigma =
                           runif(1, 0, 5), beta = runif(1, -5, 5), mean.p = runif(1, 0, 1))}

# Parameters monitored
parameters <- c("mean.phi", "mean.p", "phi.est", "sigma2", "beta")

# MCMC settings
ni <- 20000
nt <- 6
nb <- 10000
nc <- 3

# Call WinBUGS from R (BRT 12 min)
cjs.cov <- bugs(bugs.data, inits, parameters, "cjs-cov-raneff.bug",
                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE,
                bugs.directory = bugs.dir, working.directory = getwd())

# Summarize posteriors
print(cjs.cov, digits = 3)

# Produce graph
par(mfrow = c(1, 2), las = 1)
hist(cjs.cov$sims.list$beta, nclass = 25, col = "gray", main = "",
     xlab = expression(beta), ylab = "Frequency")
abline(v = -0.3, col = "red", lwd = 2)
hist(cjs.cov$sims.list$sigma2, nclass = 50, col = "gray", main = "",
     xlab = expression(sigma^2), ylab = "Frequency", xlim=c(0, 3))
abline(v = 0.2, col = "red", lwd = 2)


# Fixed Group Effects -----------------------------------------------------

# Define parameter values
n.occasions <- 12 # Number of capture occasions
marked <- rep(30, n.occasions-1) # Annual number of newly marked
individuals
phi.f <- rep(0.65, n.occasions-1) # Survival of females
p.f <- rep(0.6, n.occasions-1) # Recapture prob. of females
phi.m <- rep(0.8, n.occasions-1) # Survival of males
p.m <- rep(0.3, n.occasions-1) # Recapture prob. of males

# Define matrices with survival and recapture probabilities
PHI.F <- matrix(phi.f, ncol = n.occasions-1, nrow = sum(marked))
P.F <- matrix(p.f, ncol = n.occasions-1, nrow = sum(marked))
PHI.M <- matrix(phi.m, ncol = n.occasions-1, nrow = sum(marked))
P.M <- matrix(p.m, ncol = n.occasions-1, nrow = sum(marked))

# Simulate capture-histories
CH.F <- simul.cjs(PHI.F, P.F, marked)
CH.M <- simul.cjs(PHI.M, P.M, marked)

# Merge capture-histories by row
CH <- rbind(CH.F, CH.M)

# Create group variable
group <- c(rep(1, dim(CH.F)[1]), rep(2, dim(CH.M)[1]))

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Specify model in BUGS language
sink("cjs-group.bug")
cat("
model {

# Priors and constraints
for (i in 1:nind){
for (t in f[i]:(n.occasions−1)){
phi[i,t] <- phi.g[group[i]]
p[i,t] <- p.g[group[i]]
} #t
} #i

for (u in 1:g){
phi.g[u] ~ dunif(0, 1) # Priors for group-specific
survival
p.g[u] ~ dunif(0, 1) # Priors for group-specific
recapture
}

# Likelihood
for (i in 1:nind){

# Define latent state at first capture
z[i,f[i]] <- 1
for (t in (f[i]+1):n.occasions){

# State process
z[i,t] ~ dbern(mu1[i,t])
mu1[i,t] <- phi[i,t−1] * z[i,t−1]

# Observation process
y[i,t] ~ dbern(mu2[i,t])
mu2[i,t] <- p[i,t−1] * z[i,t]
} #t
} #i
}
",fill = TRUE)
sink()

# Bundle data
bugs.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2],
                  z = known.state.cjs(CH), g = length(unique(group)), group = group)

# Initial values
inits <- function(){list(z = cjs.init.z(CH, f), phi.g = runif(length
                                                              (unique(group)), 0, 1), p.g = runif(length(unique(group)), 0, 1))}

# Parameters monitored
parameters <- c("phi.g", "p.g")

# MCMC settings
ni <- 5000
nt <- 3
nb <- 2000
nc <- 3

# Call WinBUGS from R (BRT 2 min)
cjs.group <- bugs(bugs.data, inits, parameters, "cjs-group.bug",
                  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE,
                  bugs.directory = bugs.dir, working.directory = getwd())

# Summarize posteriors
print(cjs.group, digits = 3)


# Random Group Effects ----------------------------------------------------

# Priors and constraints
for (i in 1:nind){
  for (t in f[i]:(n.occasions-1)){
    logit(phi[i,t]) <- beta[group[i]]
    p[i,t] <- mean.p
  } #t
} #i

for (u in 1:g){
  beta[u] ~ dnorm(mean.beta, tau)
  phi.g[u] <- 1 / (1+exp(-beta[u])) # Back-transformed group-specific survival
}

mean.beta ~ dnorm(0, 0.001) # Prior for logit of mean survival
mean.phi <- 1 / (1+exp(-mean.beta)) # Back-transformed mean survival
sigma ~ dunif(0, 10) # Prior for sd of logit of survival variability
tau <- pow(sigma, -2)
mean.p ~ dunif(0, 1) # Prior for mean recapture

# Define parameter values
n.occasions <- 20 # Number of capture occasions
marked <- rep(30, n.occasions-1) # Annual number of newly marked individuals
mean.phi <- 0.65
p <- rep(0.4, n.occasions-1)
v.ind <- 0.5

# Draw annual survival probabilities
logit.phi <- rnorm(sum(marked), qlogis(mean.phi), v.ind^0.5)
phi <- plogis(logit.phi)

# Define matrices with survival and recapture probabilities
PHI <- matrix(phi, ncol = n.occasions-1, nrow = sum(marked),
              byrow = FALSE)
P <- matrix(p, ncol = n.occasions-1, nrow = sum(marked))

# Simulate capture-histories
CH <- simul.cjs(PHI, P, marked)

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Specify model in BUGS language
sink("cjs-ind-raneff.bug")
cat("
model {

# Priors and constraints
for (i in 1:nind){
for (t in f[i]:(n.occasions−1)){
logit(phi[i,t]) <- mu + epsilon[i]
p[i,t] <- mean.p
} #t
} #i

for (i in 1:nind){
epsilon[i] ~ dnorm(0, tau)
}

mean.phi ~ dunif(0, 1) # Prior for mean survival
mu <- log(mean.phi / (1-mean.phi)) # Logit transformation
sigma ~ dunif(0, 5) # Prior for standard deviation
tau <- pow(sigma, −2)
sigma2 <- pow(sigma, 2)
mean.p ~ dunif(0, 1) # Prior for mean recapture

# Likelihood
for (i in 1:nind){

# Define latent state at first capture
z[i,f[i]] <- 1
for (t in (f[i]+1):n.occasions){

# State process
z[i,t] ~ dbern(mu1[i,t])
mu1[i,t] <- phi[i,t-1] * z[i,t-1]

# Observation process
y[i,t] ~ dbern(mu2[i,t])
mu2[i,t] <- p[i,t−1] * z[i,t]
} #t
} #i
}
",fill = TRUE)
sink()

# Bundle data
bugs.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2],
                  z = known.state.cjs(CH))

# Initial values
inits <- function(){list(z = cjs.init.z(CH, f), mean.phi = runif(1, 0, 1),
                         mean.p = runif(1, 0, 1), sigma = runif(1, 0, 2))}

# Parameters monitored
parameters <- c("mean.phi", "mean.p", "sigma2")

# MCMC settings
ni <- 50000
nt <- 6
nb <- 20000
nc <- 3

# Call WinBUGS from R (BRT 73 min)
cjs.ind <- bugs(bugs.data, inits, parameters, "cjs-ind-raneff.bug",
                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE,
                bugs.directory = bugs.dir, working.directory = getwd())

# Summarize posteriors
print(cjs.ind, digits = 3)

# Produce graph
par(mfrow = c(1, 2), las = 1)
hist(cjs.ind$sims.list$mean.phi, nclass = 25, col = "gray", main = "",
     xlab = expression(bar(phi)), ylab = "Frequency")
abline(v = mean.phi, col = "red", lwd = 2)
hist(cjs.ind$sims.list$sigma2, nclass = 15, col = "gray", main = "",
     xlab = expression(sigma^2), ylab = "Frequency", xlim = c(0, 3))
abline(v = v.ind, col = "red", lwd = 2)

## Estimate survival as a function of an individual covariate x
# Priors and constraints
for (i in 1:nind){
  for (t in f[i]:(n.occasions-1)){
    logit(phi[i,t]) <- mu + beta*x[i] + epsilon[i]
    p[i,t] <- mean.p
  } #t
} #i

for (i in 1:nind){
  epsilon[i] ~ dnorm(0, tau)
}

mean.phi ~ dunif(0, 1) # Prior for mean survival
mu <- log(mean.phi / (1-mean.phi)) # Logit transformation
beta ~ dnorm(0, 0.001) # Prior for covariate slope
sigma ~ dunif(0, 5) # Prior for standard deviation
tau <- pow(sigma, -2)
sigma2 <- pow(sigma, 2)
mean.p ~ dunif(0, 1) # Prior for mean recapture


# Fixed Group and Time effects --------------------------------------------

# Define parameter values
n.occasions <- 12 # Number of capture occasions
marked <- rep(50, n.occasions-1) # Annual number of newly marked individuals
phi.f <- c(0.6, 0.5, 0.55, 0.6, 0.5, 0.4, 0.6, 0.5, 0.55, 0.6, 0.7)
p.f <- rep(0.6, n.occasions-1)
diff <- 0.5 # Difference between male and female survival on logit scale
phi.m <- plogis(qlogis(phi.f) + diff)
p.m <- rep(0.3, n.occasions-1)

# Define matrices with survival and recapture probabilities
PHI.F <- matrix(rep(phi.f, sum(marked)), ncol = n.occasions-1,
                nrow = sum(marked), byrow = TRUE)
P.F <- matrix(rep(p.f, sum(marked)), ncol = n.occasions-1,
              nrow = sum(marked), byrow = TRUE)
PHI.M <- matrix(rep(phi.m, sum(marked)), ncol = n.occasions-1,
                nrow = sum(marked), byrow = TRUE)
P.M <- matrix(rep(p.m, sum(marked)), ncol = n.occasions-1,
              nrow = sum(marked), byrow = TRUE)

# Simulate capture-histories
CH.F <- simul.cjs(PHI.F, P.F, marked)
CH.M <- simul.cjs(PHI.M, P.M, marked)

# Merge capture-histories
CH <- rbind(CH.F, CH.M)

# Create group variable
group <- c(rep(1, dim(CH.F)[1]), rep(2, dim(CH.M)[1]))

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Specify model in BUGS language
sink("cjs-add.bug")
cat("
model {

# Priors and constraints
for (i in 1:nind){
for (t in f[i]:(n.occasions−1)){
logit(phi[i,t]) <- beta[group[i]] + gamma[t]
p[i,t] <- p.g[group[i]]
} #t
} #i

# for survival parameters
for (t in 1:(n.occasions−1)){
gamma[t] ~ dnorm(0, 0.01)I(−10, 10) # Priors for time effects
phi.g1[t] <- 1 / (1+exp(−gamma[t])) # Back-transformed survival of males
phi.g2[t] <- 1 / (1+exp(−gamma[t]-beta[2])) # Back-transformed survival of females
}

beta[1] <- 0 # Corner constraint
beta[2] ~ dnorm(0, 0.01)I(−10, 10) # Prior for difference in male and female survival

# for recapture parameters
for (u in 1:g){
p.g[u] ~ dunif(0, 1) # Priors for group-spec. recapture
}

# Likelihood
for (i in 1:nind){

# Define latent state at first capture
z[i,f[i]] <- 1
for (t in (f[i]+1):n.occasions){

# State process
z[i,t] ~ dbern(mu1[i,t])
mu1[i,t] <- phi[i,t−1] * z[i,t−1]

# Observation process
y[i,t] ~ dbern(mu2[i,t])
mu2[i,t] <- p[i,t−1] * z[i,t]
} #t
} #i
}
",fill = TRUE)
sink()

# Bundle data
bugs.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2],
                  z = known.state.cjs(CH), g = length(unique(group)), group = group)

# Initial values
inits <- function(){list(z = cjs.init.z(CH, f), gamma =
                           rnorm(n.occasions-1), beta = c(NA, rnorm(1)), p.g = runif(length
                                                                                     (unique(group)), 0, 1))}

# Parameters monitored
parameters <- c("phi.g1", "phi.g2", "p.g", "beta")

# MCMC settings
ni <- 5000
nt <- 3
nb <- 2000
nc <- 3

# Call WinBUGS from R (BRT 7 min)
cjs.add <- bugs(bugs.data, inits, parameters, "cjs-add.bug",
                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE,
                bugs.directory = bugs.dir, working.directory = getwd())

# Summarize posteriors
print(cjs.add, digits = 3)

# Figure of male and female survival
lower.f <- upper.f <- lower.m <- upper.m <- numeric()
for (t in 1:(n.occasions-1)){
  lower.f[t] <- quantile(cjs.add$sims.list$phi.g1[,t], 0.025)
  upper.f[t] <- quantile(cjs.add$sims.list$phi.g1[,t], 0.975)
  lower.m[t] <- quantile(cjs.add$sims.list$phi.g2[,t], 0.025)
  upper.m[t] <- quantile(cjs.add$sims.list$phi.g2[,t], 0.975)
}
plot(x=(1:(n.occasions-1))-0.1, y = cjs.add$mean$phi.g1, type = "b",
     pch = 16, ylim = c(0.2, 1), ylab = "Survival probability",
     xlab = "Year", bty = "n", cex = 1.5, axes = FALSE)
axis(1, at = 1:11, labels = rep(NA,11), tcl = -0.25)
axis(1, at = seq(2,10,2), labels = c("2","4","6","8","10"))
axis(2, at = seq(0.2, 1, 0.1), labels = c("0.2", NA, "0.4", NA, "0.6", NA,
                                          "0.8", NA, "1.0"), las = 1)
segments((1:(n.occasions-1))-0.1, lower.f, (1:(n.occasions-1))-0.1,
         upper.f)
points(x = (1:(n.occasions-1))+0.1, y = cjs.add$mean$phi.g2,
       type = "b", pch = 1, lty = 2, cex = 1.5)
segments((1:(n.occasions-1))+0.1, lower.m, (1:(n.occasions-1))+0.1,
         upper.m)

## Fit model with interaction of sex and time
# Priors and constraints
for (i in 1:nind){
  for (t in f[i]:(n.occasions-1)){
    phi[i,t] <- eta.phi[group[i],t]
    p[i,t] <- p.g[group[i]]
  } #t
} #i
# for survival parameters
for (u in 1:g){
  for (t in 1:(n.occasions-1)){
    eta.phi[u,t] ~ dunif(0, 1) # Prior for time and group-spec. survival
  } #t
} #g

# for recapture parameters
for (u in 1:g){
  p.g[u] ~ dunif(0, 1) # Priors for group-spec. recapture
}



# Fixed Group and Random Time effects -------------------------------------

# Priors and constraints
for (i in 1:nind){
  for (t in f[i]:(n.occasions-1)){
    logit(phi[i,t]) <- eta.phi[group[i],t]
    p[i,t] <- p.g[group[i]]
  } #t
} #i

# for survival parameters
for (u in 1:g){
  for (t in 1:(n.occasions-1)){
    eta.phi[u,t] <- mu.phi[u] + epsilon[u,t]
    epsilon[u,t] ~ dnorm(0, tau[u])
  } #t
  
  mean.phi[u] ~ dunif(0, 1) # Priors on mean group-spec.
  survival
  mu.phi[u] <- log(mean.phi[u] / (1-mean.phi[u]))
  sigma[u] ~ dunif(0, 10) # Priors for group-spec. sd
  tau[u] <- pow(sigma[u], -2)
  sigma2[u] <- pow(sigma[u], 2)
} #g

# for recapture parameters
for (u in 1:g){
  p.g[u] ~ dunif(0,1) # Priors for group-spec.recapture
}

# Specify model in BUGS language
sink("cjs-temp-corr.bug")
cat("
model {

# Priors and constraints
for (i in 1:nind){
for (t in f[i]:(n.occasions−1)){
logit(phi[i,t]) <- eta.phi[t,group[i]]
p[i,t] <- p.g[group[i]]
} #t
} #i

# for survival parameters
for (t in 1:(n.occasions−1)){
eta.phi[t,1:g] ~ dmnorm(mu.phi[], Omega[,])
} #t

for (u in 1:g){
mean.phi[u] ~ dunif(0, 1) # Priors on mean group-spec. survival
mu.phi[u] <- log(mean.phi[u] / (1−mean.phi[u]))
} #g

Omega[1:g, 1:g] ~ dwish(R[,], df) # Priors for variance-covariance
matrix
Sigma[1:g, 1:g] <- inverse(Omega[,])

# for recapture parameters
for (u in 1:g){
p.g[u] ~ dunif(0, 1) # Priors for group-spec. recapture
}

# Likelihood
for (i in 1:nind){

# Define latent state at first capture
z[i,f[i]] <- 1
for (t in (f[i]+1):n.occasions){

# State process
z[i,t] ~ dbern(mu1[i,t])
mu1[i,t] <- phi[i,t−1] * z[i,t−1]

# Observation process
y[i,t] ~ dbern(mu2[i,t])
mu2[i,t] <- p[i,t−1] * z[i,t]
} #t
} #i
}
",fill = TRUE)
sink()

# Bundle data
bugs.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2],
                  z = known.state.cjs(CH), g = length(unique(group)), group = group,
                  R = matrix(c(5, 0, 0, 1), ncol = 2), df = 3)

# Initial values
inits <- function(){list(z = cjs.init.z(CH, f), p.g = runif(length
                                                            (unique(group)), 0, 1), Omega = matrix(c(1, 0, 0, 1), ncol = 2))}

# Parameters monitored
parameters <- c("eta.phi", "p.g", "Sigma", "mean.phi")

# MCMC settings
ni <- 5000
nt <- 3
nb <- 2000
nc <- 3

# Call WinBUGS from R (BRT 5 min)
cjs.corr <- bugs(bugs.data, inits, parameters, "cjs-temp-corr.bug",
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE,
                 bugs.directory = bugs.dir, working.directory = getwd())

# Summarize posteriors
print(cjs.corr, digits = 3)

corr.coef <- cjs.corr$sims.list$Sigma[,1,2] / sqrt(cjs.corr$sims.list$Sigma[,1,1] * cjs.corr$sims.list$Sigma[,2,2])
