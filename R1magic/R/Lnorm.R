#'
#' L-p norm of a given complex vector
#'
#'@author Mehmet Suzen
#'@param X, a complex vector, can be real too.
#'@param p, norm value
#'@return L-p norm of the complex vector
#'
Lnorm <- function(X, p) {
  LX <- Mod(X)
  if(p>0 & p <Inf) {
    LX <- sum(LX^p)^(1.0/p)
  }
  if(p == Inf) return(max(LX))
  if(p == 0) return(length(which(LX==0)))
  return(LX)
}
