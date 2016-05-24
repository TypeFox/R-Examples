
#'  Compute Guttman's Lambda 5 Coefficient
#'
#' @param x Can be either a data matrix or a covariance matrix.
#' @param missing how to handle missing values.
#' @param standardize Results are standardized by using the correlation matrix instead of the covariance matrix for computation.
#'
#' @references
#' Guttman L (1945). "A Basis for Analyzing Test-Retest Reliability." Psychometrika, 10, 255-282.
#' @author Tyler Hunt \email{tyler@@psychoanalytix.com}
#' @examples
#' lambda5(Rosenberg)
#' @export
lambda5 <- function(x, missing="complete", standardize=FALSE){
  
  n <- dim(x)[1]
  p <- dim(x)[2]
    
  sigma <- impute.cov(x, missing)
  
  if(standardize==TRUE){
    sigma <- cov2cor(sigma)
  }
    
  sigma0 <- sigma
  diag(sigma0) <- 0
  items <- nrow(sigma)
	
  t1  <- matrix(rep(1, items), nrow=1)
  t1t <- t(t1)
	L1  <- 1-sum(diag(sigma))/sum(sigma)
	gamma2max <- max(colSums(sigma0^2))
	lambda5 <- L1 + 2*sqrt(gamma2max)/sum(sigma)

    structure(c(lambda5=lambda5), class = "lambda5")
}
