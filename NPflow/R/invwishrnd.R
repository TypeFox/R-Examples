#' Sample from a inverse-Wishart distribution
#' 
#' For internal use only.
#' 
#'@param n degrees of freedom
#'
#'@param lambda scale parameter
#' 
#'@keywords internal
#'
#'@export
#'
#'
invwishrnd <- function(n,lambda){
  p<-ncol(lambda)
  S<-try(solve(wishrnd(n = n, Sigma = solve(lambda))),silent=TRUE)
  if(inherits(S, "try-error")){
    S=solve(wishrnd(n = n, Sigma = solve((lambda+diag(p)))))
  }
  return(S)
}
