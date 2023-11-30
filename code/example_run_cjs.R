
# setup -------------------------------------------------------------------

pacman::p_load(tidyverse,
               lme4)

source("code/function.R")

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
## n0: number of tagging per occasion
## m_occ: number of occasion

obj <- cjs_data()