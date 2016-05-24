NULL
#'
#' This function finds the bivariate joint probability or the binary correlation from the corresponding Gaussian correlation \code{x}
#'
#' @param x value of expected correlation between the corresponding Gaussian-distributed variables 
#' @param p0_v1,p0_v2 probability of no precipitation occurences for the v1 and v2 time series respectively. See \code{Notes}.
#' @param correlation logical numeric value. Default is \code{FALSE}. If \code{TRUE} the function returns the binary correlation like eq. 6 of Mhanna, et al.,2011.
# @param ... 
#' 
#' @author Emanuele Cordano
#' 
#' @return probability of no precipitation occurence in both v1 and v2 simultaneously. It is a matrix if \code{x} is a matrix.
#' 
#' @note This function makes use of normal copula. A graphical introduction to this function (with its inverse) makes is present in the following URL references:  \url{http://onlinelibrary.wiley.com/doi/10.1002/joc.2305/abstract} 
#'   and \url{http://www.sciencedirect.com/science/article/pii/S0022169498001863} (See fig. 1 and par. 3.2) 
#' If the argument \code{p0_v2}, the two marginal probabily values must be given as a vector through the  argument \code{p0_v1}: \code{p0_v1=c(p0_v1,p0_v2)} . 
#' In case \code{x} is a correlation/covariance matrix the marginal probabilities are given as a vector through the argument \code{p0_v1}.
#'
#' @references
#' D.S. Wilks (1998), Multisite Generalization of a Daily Stochastic Precipitation Generation Model, Journal of Hydrology, Volume 210, Issues 1-4, September 1998, Pages 178-191,
#' \url{http://www.sciencedirect.com/science/article/pii/S0022169498001863}
#' 
#' Muamaraldin Mhanna and Willy Bauwens (2011) A Stochastic Space-Time Model for the Generation of Daily Rainfall in the Gaza Strip, International Journal of Climatology, Volume 32, Issue 7, pages 1098-1112,
#' \url{http://dx.doi.org/10.1002/joc.2305}
#' 
#'  
#'  
#' 
#' @seealso \code{\link{normalCopula}},\code{\link{pcopula}}

#' @import copula
#' @export
#' @examples 

#' rho <- 0.4
#' p00 <- omega(x=rho,p0_v1=0.5,p0_v2=0.5)
#' cor00 <- omega(x=rho,p0_v1=0.5,p0_v2=0.5,correlation=TRUE)

## str

omega <- function(x=0.5,p0_v1=0.5,p0_v2=NA,correlation=FALSE) {
	
	out <- NA 
	if (is.na(p0_v2)) {
		if (length(p0_v1)>1) { 
			p0_v2 <- p0_v1[2]
		} else {
			p0_v2 <- p0_v1
		}
		
	}
	if (is.matrix(x)) {
		NR <- nrow(x)
		NC <- ncol(x)
		out <- array(NA,c(NR,NC))
		for (r in 1:NR) {
			for (c in 1:NC) {
				
				out[r,c] <- omega(x=x[r,c],p0_v1=p0_v1[r],p0_v2=p0_v1[c],correlation=correlation)
				
			}
		}
	
		return(out)
		
	}
	
	if (length(x)==1) {
	
		nc <- normalCopula(x,dim=2) ## t is a considered a 2d Normal Copula!!! #### 
		out <- pCopula(copula=nc,u=c(p0_v1,p0_v2))
	} else if (length(x)>1) {
		
		out <- as.numeric(lapply(x,omega,p0_v1=p0_v1,p0_v2=p0_v2))
		
	}
	
	if (correlation) {
		
		out <- (out-p0_v1*p0_v2)/(p0_v1*p0_v2*(1-p0_v1)*(1-p0_v2))^0.5
		
	}
	return(out)
}