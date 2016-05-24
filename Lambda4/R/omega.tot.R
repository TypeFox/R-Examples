#' Compute McDonald's Omega Total
#' 
#' @description  McDonald proposed Omega Total as a method for estimating reliabilty for a test with multiple factors.
#'
#' @param x Can be either a data matrix or a covariance matrix
#' @param missing how to handle missing values. \eqn{mi}.
#' @param factors The number of latent factors.
#' 
#' @return
#' \item{omega.tot}{Omega total reliability estimate.}
#' 
#' @author Tyler Hunt \email{tyler@@psychoanalytix.com}
#' @references
#' McDonald, R. P. (1999). Test Theory: Aunified Treatment. Psychology Press.
#' @examples
#' omega.tot(Rosenberg, factors=1)
#' 
#' @export
omega.tot<-function(x, factors=1, missing="complete"){

	stopifnot(require(GPArotation))

	n <- dim(x)[1]
	p <- dim(x)[2]
	sigma <- impute.cov(x, missing)
  
  sigma<-cov2cor(sigma)
  items<-nrow(sigma)
	
    ff<-factanal(covmat=sigma, factors=factors)

    omega.tot<- 1-sum(ff$uniquenesses)/sum(sigma)
    result<-c(omega.tot=omega.tot)
    class(result)<-c("omega.tot")
    return(result)
}

