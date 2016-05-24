#' Compute Guttman's Lambda 1 Coefficient
#'
#' @param x an object that you can compute the covariance of
#' @param missing how to handle missing values.
#' @param standardize Results are standardized by using the correlation matrix instead of the covariance matrix for computation.
#'
#' @references
#' Guttman L (1945). "A Basis for Analyzing Test-Retest Reliability." Psychometrika, 10, 255-282.
#' @author Tyler Hunt \email{tyler@@psychoanalytix.com}
#' @examples
#' lambda1(Rosenberg)
#' @export
lambda1<-function(x, missing="complete", standardize=FALSE){
  
  n <- dim(x)[1]
  p <- dim(x)[2]
  
  sigma <- impute.cov(x, missing)
  
  if(standardize==TRUE){
    sigma <- cov2cor(sigma)
  }
	
	t1  <- matrix(rep(1, nrow(sigma)), nrow=1)
	t1t <- t(t1)

	lambda1 <- 1-sum(diag(sigma))/sum(sigma)
    
    structure(c(lambda1=lambda1), class = "lambda1")
}
