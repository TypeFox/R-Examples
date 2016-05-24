.tsd_fit <- function(par, fixed.parameters=NULL, males, N, temperatures, equation) {

  par <- c(par, fixed.parameters)
  p <- .modelTSD(par, temperatures, equation)
  p <- ifelse(p==0, 1E-9, p)
  p <- ifelse(p==1, 1-1E-9, p)

    if (any(is.infinite(p))) {return(Inf)} else {
   return(-sum(dbinom(males, N, p, log = TRUE)))
  }
  
}
