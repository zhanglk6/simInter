SRS <- function(n, pi) {
  ## simple randomization
  sample(c(0,1), n, T, c(1-pi, pi))
}

one_block <- function(pi, bsize){
  ## sample a block 
  num_of_1 = floor(bsize*pi)
  num_of_0 = bsize - num_of_1
  sample(c(rep(1, num_of_1), rep(0, num_of_0)), bsize, F)
}

SBR <- function(strata, pi, bsize = 6) {
  ## stratified block randomization
  ## strata: a vector of strata of size n
  strt_num = max(strata)
  n = length(strata)
  counting  = numeric(length = strt_num)
  T_STR = numeric(length = n)
  BG = sapply(1:strt_num, function(n) one_block(pi,bsize))
  
  for (i in 1:n) {
    r = strata[i]
    counting[r] = counting[r] + 1
    pos = counting[r] %% bsize
    if(pos == 0){
      T_STR[i] = BG[bsize, r] 
      BG[,r] = one_block(pi,bsize)
    }
    else{
      T_STR[i] = BG[pos, r]
    }
  }
  return(T_STR)
}


SBCD <- function(strata, pi, lambda = 0.75) {
  ## stratified Efron's biased coin design
  strt_num = max(strata)
  n = length(strata)
  D  = numeric(length = strt_num)
  T_STR = numeric(length = n)
  
  brid = c(1.0/pi, -1.0/(1.0-pi))
  
  u = runif(n)
  
  for (i in 1:n) {
    r = strata[i]
    argval = D[r]
    if(argval > -0.000001 && argval < 0.000001){
     # T_one = u[i] < 0.5 change in 2024.09.03
      T_one = u[i] < pi
    } else if (argval > 0.000001) {
      T_one = u[i] < 1 - lambda
    } else {
      T_one = u[i] < lambda
    }
    D[r] = D[r] + brid[2 - T_one]
    T_STR[i] = T_one
  }
  return(T_STR)
}

PS <- function(profiles, weight, pi, lambda = 0.75) {
  ## stratified Efron's biased coin design
  n_cov = ncol(profiles)
  n = nrow(profiles)
  cov_level_num = c(0,apply(profiles, 2, max, na.rm = TRUE) + 1) ## first position for reservation
  D  = numeric(length = sum(cov_level_num))
  r = numeric(n_cov)
  T_STR = numeric(length = n)
  
  brid = c(1.0/pi, -1.0/(1.0-pi))
  
  u = runif(n)
  
  for (i in 1:n) {
    profile = profiles[i,]
    r = profile + 1 + cov_level_num[1:n_cov]
    argval = sum(weight * D[r])
    if(argval > -0.000001 && argval < 0.000001){
      T_one = u[i] < pi
    } else if (argval > 0.000001) {
      T_one = u[i] < 1 - lambda
    } else {
      T_one = u[i] < lambda
    }
    D[r] = D[r] + brid[2 - T_one]
    T_STR[i] = T_one
  }
  return(T_STR)
}

AdaptFUN <- function(d, ns, pi) {
  if(ns == 0 || d == 0){
    return(pi)
  } else if(d > 0){
    return(pi*(1.0 - d/ns))
  } else{
    return((pi-1.0)*d/ns + pi)
  }
}

WEI <- function(strata, pi) {
  ## WEI's urn design
  strt_num = max(strata)
  n = length(strata)
  D  = numeric(length = strt_num)
  nk  = numeric(length = strt_num)
  T_STR = numeric(length = n)
  
  brid = c(1.0/pi, -1.0/(1.0-pi))
  
  u = runif(n)
  
  for (i in 1:n) {
    r = strata[i]
    p = AdaptFUN(D[r],nk[r],pi)
    T_one = u[i] < p
    D[r] = D[r] + brid[2 - T_one]
    nk[r] = nk[r] + abs(brid[2 - T_one])
    T_STR[i] = T_one
  }
  return(T_STR)
}
