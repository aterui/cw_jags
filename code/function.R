
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



# mrcheck function --------------------------------------------------------

mrcheck <- function(data,
                    xi = 0.3,
                    rm_multi = TRUE,
                    cnm = c("occasion",
                            "date",
                            "time",
                            "tag_id",
                            "site",
                            "species",
                            "section",
                            "length",
                            "weight",
                            "recap",
                            "mortality",
                            "fin_clip",
                            "fin_recap",
                            "glued")) {
  
  pacman::p_load(tidyverse,
                 foreach)
  
  colnames(data) <- str_to_lower(colnames(data))
  
  ## check column input  
  if (any(colnames(data) == "tag_id2"))
    stop("tag_id2 exists; replace tag_id with tag_id2")
  
  if(sum(colnames(data) %in% cnm) != length(cnm))
    stop(paste("Missing column(s)",
               sQuote(paste(setdiff(cnm, colnames(data)),
                            collapse = ", "))
    )
    )
  
  ## drop rows with no tag_id
  if (any(is.na(data$tag_id))) {
    data <- drop_na(data, tag_id)
    message("Rows with NA tag_id were dropped")
  }
  
  ## check species consistency
  df_n <- data %>% 
    group_by(tag_id) %>% 
    summarize(n = n_distinct(species)) %>% 
    filter(n > 1)
  
  if (nrow(df_n) > 0) {
    message("Species ID check - ERROR: One or more tags have multiple species ID")
    return(df_n)
  } else {
    message("Species ID check - PASS: Each tag has unique species ID")
  }
  
  ## size - weight relationship
  sp_i <- with(data, unique(species) %>%
                 sort() %>%
                 na.omit())
  
  list_df <- foreach(i = seq_len(length(sp_i))) %do% {
    
    df_s <- filter(data, species == sp_i[i])
    
    m <- MASS::rlm(log(weight) ~ log(length) + factor(occasion),
                   df_s)
    
    return(mutate(df_s, rlm_weight = m$w))
  }
  
  df_w <- do.call(bind_rows, list_df)  
  
  if (any(with(df_w, rlm_weight) < xi)) {
    message("One or more data points have suspicious length or weight entries; check column 'rlm_weight'")
    df_q <- filter(df_w, rlm_weight < xi)
    return(df_q)
  }
  
  ## multiple captures of the same individuals within an occasion
  if (rm_multi) {
    df_multi <- df_w %>% 
      group_by(tag_id, occasion) %>% 
      summarize(n_obs = n()) %>% 
      ungroup() %>% 
      filter(n_obs > 1) %>% 
      arrange(occasion) %>% 
      data.table::data.table()
    
    if (nrow(df_multi) > 0) {
      message("Multiple (re)captures within an occasion exist; only last captures were retained")
      df_out <- df_w %>% 
        mutate(datetime = paste(date, time),
               datetime = as.POSIXct(datetime,
                                     format = "%m/%d/%Y %H:%M:%OS")) %>% 
        relocate(datetime) %>% 
        group_by(tag_id, occasion) %>% 
        slice(which.max(datetime)) %>% 
        ungroup()
      
      attr(df_out, 'multi_recap') <- df_multi
    } else {
      df_out <- df_w
    }
    
    return(df_out)
  } else {
    message("No changes were made to the data")
    return(data)
  }
  
}



# get first ---------------------------------------------------------------

getf <- function(x) min(which(x == 1))
