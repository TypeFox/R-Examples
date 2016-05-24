BiCopTau2Par <- function(family, tau) {
  
  if (!(family %in% c(0, 1, 3, 4, 5, 6, 13, 14, 16, 23, 24, 26, 33, 34, 36))) 
    stop("Copula family not implemented.")
  
  if (family == 0) {
    par <- 0
  } else if (family == 1) {
    par <- sin(pi * tau/2)
  } else if (family %in% c(3, 13)) {
    if (tau <= 0) 
      stop("Clayton copula cannot be used for tau<=0.")
    par <- 2 * tau/(1 - tau)
  } else if (family %in% c(4, 14)) {
    if (tau < 0) 
      stop("Gumbel copula cannot be used for tau<0.")
    par <- 1/(1 - tau)
  } else if (family == 5) {
    if (tau == 0) 
      stop("Frank copula cannot be used for tau=0.")
    par <- Frank.itau.JJ(tau)
  } else if (family %in% c(6, 16)) {
    if (tau <= 0) 
      stop("Joe copula cannot be used for tau<=0.")
    par <- Joe.itau.JJ(tau)
  } else if (family %in% c(23, 33)) {
    if (tau >= 0) 
      stop("Rotated Clayton copula cannot be used for tau>=0.")
    par <- 2 * tau/(1 + tau)
  } else if (family %in% c(24, 34)) {
    if (tau > 0) 
      stop("Rotated Gumbel copula cannot be used for tau>0.")
    par <- -(1/(1 + tau))
  } else if (family %in% c(26, 36)) {
    if (tau >= 0) 
      stop("Rotated Joe copula cannot be used for tau>=0.")
    par <- -Joe.itau.JJ(-tau)
  }
  
  return(par)
  
} 
