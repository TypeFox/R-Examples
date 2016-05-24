#' Compute Guttman's Lambda 2 Coefficient
#'
#' @param x Can be either a data matrix or a covariance matrix
#' @param missing how to handle missing values.
#' @param standardize Results are standardized by using the correlation matrix instead of the covariance matrix for computation.
#'
#' @references
#' Guttman L (1945). "A Basis for Analyzing Test-Retest Reliability." Psychometrika, 10, 255-282.
#' @author Tyler Hunt \email{tyler@@psychoanalytix.com}
#' @examples
#' lambda2(Rosenberg)
#' @export
lambda2<-function(x, missing="complete", standardize=FALSE){
  
  n <- dim(x)[1]
  p <- dim(x)[2]
  
  sigma <- impute.cov(x, missing)
  
  if(standardize==TRUE){
    sigma <- cov2cor(sigma)
  }

	sigma0<-sigma
	diag(sigma0)=0
	
	t1<-matrix(rep(1, n), nrow=1)
	t1t<-t(t1)
	
	gamma2<-sum(sigma0^2)
	L1<-1-sum(diag(sigma))/sum(sigma)
	
	lambda2<-L1+sqrt((n/(n-1))*gamma2)/sum(sigma)

result<-c(lambda2=lambda2)
class(result)<-c("lambda2")
return(result)
}
