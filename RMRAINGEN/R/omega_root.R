NULL
#'
#' This is the target function whose zero is searched to crete the inverse function of  \code{\link{omega}}.
#'
#' @param x value of expected correlation between the corresponding Gaussian-distributed variables 
#' @param p0_v1,p0_v2 probablity of no precipitatin occurences for the v1 and v2 time series respectively. 
#' @param p00 probability of no precipitation occurence in both v1 and v2 simultanously returned by \code{\link{omega}}
#' @param correlation numerical value. DEfault is \code{NA}.  Binary correlation retured by \code{\link{omega}}  when the argumet \code{correlation=TRUE}
#' 
#' @author Emanuele Cordano
#' 
#' @return  the value \code{p00-omega(x=x,p0_v1=p0_v1,p0_v2=p0_v2)} or \code{correlation-omega(x=x,p0_v1=p0_v1,p0_v2=p0_v2)} (if \code{correlation} is not \code{NA})
#' 
#' @note This function makes use of normal copula 
#' 
#' @seealso \code{\link{normalCopula}},\code{\link{pcopula}},\code{\link{omega}},\code{\link{omega_inv}}
#' @export
#' @examples 

#' rho <- 0.4
#' p00 <- omega(x=rho,p0_v1=0.5,p0_v2=0.5)
#' omega_root(x=rho,p0_v1=0.5,p0_v2=0.5,p00=p00)

# TO GO ON with omega ....

omega_root <- function(x=0.5,p0_v1=0.5,p0_v2=0.5,p00=p0_v1*p0_v2,correlation=NA) {
	
	
	if (is.na(correlation)) {
		out <- p00-omega(x,p0_v1=p0_v1,p0_v2=p0_v2,correlation=FALSE)
	} else { 
		out <- 	correlation-omega(x,p0_v1=p0_v1,p0_v2=p0_v2,correlation=TRUE)
	}
	return(out)
}