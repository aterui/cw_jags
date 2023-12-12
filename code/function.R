
# simulation function for cjs model ---------------------------------------
## s: survival probability
## p: detection probability
## n0: number of tagging per occasion
## m_occ: number of occasion

cjs_data <- function(s = 0.9, p = 0.5, n0 = 30, n_occ = 10) {
  
  m_z <- matrix(NA, nrow = n0 * n_occ, ncol = n_occ)  # m_z: matrix of true survival state
  m_y <- matrix(NA, nrow = n0 * n_occ, ncol = n_occ)  # m_y: matrix of observed "(re)capture" state
  
  ## initial tag condition
  m_z[1:n0, 1] <- 1
  m_y[1:n0, 1] <- 1
  
  ## tagging over time  
  for (t in 1:(n_occ - 1)) {
    ## count tagged individuals at time t
    n_tag <- sum(!is.na(m_z[, t]))
    
    ## survival process of tagged inds 
    m_z[1:n_tag, t + 1] <- rbinom(n = n_tag, 
                                  size = m_z[1:n_tag, t],
                                  prob = s)
    m_z[(n_tag + 1):(n_tag + n0), t + 1] <- 1
    
    ## observation process of imperfect detection
    phi <- rbinom(n = n_tag, size = 1, prob = p) 
    m_y[1:n_tag, t + 1] <- m_z[1:n_tag, t + 1] * phi
    m_y[(n_tag + 1):(n_tag + n0), t + 1] <- 1
  }
  
  return(list(Z = m_z,
              Y = m_y))
}
