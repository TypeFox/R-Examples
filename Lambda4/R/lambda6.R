#' Compute Guttman's Lambda 6 Coefficient
#'
#' @param x Can be either a data matrix or a covariance matrix.
#' @param missing how to handle missing values.
#' @param standardize Results are standardized by using the correlation matrix instead of the covariance matrix for computation.
#'
#' @references
#' Guttman L (1945). "A Basis for Analyzing Test-Retest Reliability." Psychometrika, 10, 255-282.
#' @author Tyler Hunt \email{tyler@@psychoanalytix.com}
#' lambda6(Rosenberg)
#' @export
lambda6<-function(x, missing="complete", standardize=FALSE){
  
  n <- dim(x)[1]
  p <- dim(x)[2]
  
  sigma <- impute.cov(x, missing)
  
  if(standardize==TRUE){
    sigma <- cov2cor(sigma)
  }
  
  smc <- 1 - 1/diag(solve(sigma))
  
  t1 <- matrix(rep(1, p), ncol=1)
  t1t <- t(t1)
  
  lambda6 <- 1 - sum(1-smc)/(sum(sigma))
  
  
  result <- c(lambda6=lambda6)
  class(result) <- c("Lambda4.pkg","lambda6")
  return(result)
}
