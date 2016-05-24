#
#
NULL
#'
#' This function is the inverse of \code{\link{omega}} function 
#'
# @param x value of expected correlation between the corresponding Gaussian-distributed variables 
#' @param p0 matrix of joint probabilities. Default is \code{NULL}, otherwise functions returns a matrix with values 
#' @param p0_v1,p0_v2 probablity of no precipitatin occurrences for the v1 and v2 time series respectively. 
#' @param p00 probability of no precipitation occurence in both v1 and v2 simultanously returned by \code{\link{omega}}
#' @param only.value logical value. If \code{TRUE} (Default) the only Gaussian correletion (\code{x} input variable of \code{\link{omega}}) is returned, 
#' otherwise  the complete output of \code{\link{uniroot}} is returned.
#' @param correlation numerical value. DEfault is \code{NA}.  Binary correlation retured by \code{\link{omega}}  when the argumet \code{correlation=TRUE} (see \code{\link{omega_root}})
#' @param interval see \code{interval} option of \code{\link{uniroot}}. Default is \code{c(-1,1)}. 
#' @param tolerance tolerance (numeric) parameter used for comparisons with the extreme value of marginal probabilities. Default is 0.001. 
#' @param nearPD logical. If \code{TRUE} (Default) a positive-definite correlation matrix is returned by applying \code{\link{nearPD}} in case \code{p0} is a matrix and not \code{NULL}.
#' @param force.independence logical value. Default is \code{TRUE}. If it is \code{TRUE}, no negative corelation is considered and negative values of correletion are forced to be 0 (independence).
#' @param ... further arguments for \code{\link{uniroot}}
#' 
#' @author Emanuele Cordano
#' 
#' @return  value of expected correlation between the corresponding Gaussian-distributed variables (see \code{x} input argument of \code{\link{omega}}.
#' 
#' @note This function finds the zero of the  \code{\link{omega_root}} function by calling \code{\link{uniroot}}. 
#' If the argument \code{p0} is not \code{NULL} and is a matrix of joint probabilities, the function returns a correlation matrix by using the elements of \code{p0} ass joint probabilities for each couple and \code{p0_v1} as a vector of marginal probability of each occurence/no-occurence
#' (In this case if the length of \code{p0_v1} does not correspond to the number of columns of \code{p0}, the marginal probabilities are taken from the diagonal of \code{p0}).
#' See the R code for major details. 
#' 
#' 
#' @import Matrix
#' @seealso \code{\link{normalCopula}},\code{\link{pcopula}},\code{\link{omega}}(and reference URLs therein)
#' @export
#' @examples 

#' 
#' x <- omega_inv(p0_v1=0.5,p0_v2=0.5,p00=1.1*0.5*0.5)
#' omega(x,p0_v1=0.5,p0_v2=0.5)






omega_inv <- function(p0=NULL,p0_v1=0.5,p0_v2=p0_v1,p00=p0_v1*p0_v2,correlation=NA,only.value=TRUE,interval=c(-1,1),tolerance=0.001,nearPD=TRUE,force.independence=TRUE,...){

   ### Insertion of force.inpendence ec 20130310	
#### force.independence=0.05	
#	if (is.null(p0) & (force.independence>0)) {
#		
#		pind <- p0_v1*p0_v2
#		pint <- abs(p00-pind)
#		p00[pint<=force.independence] <- pind
#		
#	}
	
	### END Insertion of force.inpendence ec 20130310	
	
	
	if (!is.null(p0)) {
		out <- NULL
		message("Hmm... p0 - first argument - must be a matrix of probabilities!!!" )
	
		out <- array(NA,c(nrow(p0),ncol(p0)))
		if (length(p0_v1)!=nrow(p0)) p0_v1 <- diag(p0)
		
		for (r in 1:nrow(p0)) {
			for (c in 1:ncol(p0)) {
			
				
				
				
				
				out[r,c] <- omega_inv(p0_v1=p0_v1[r],p0_v2=p0_v1[c],p00=p0[r,c],correlation=NA,only.value=TRUE,interval=interval,tolerance=tolerance,force.independence=TRUE,...)
		
			
			
			}
			
		}
		
		if (nearPD) {
			
			out[(out<0) & (force.independence)] <- 0
			
			out <- (nearPD(out,keepDiag=TRUE,maxit=100000)) ## mod by 
			
			out <- as.matrix(out$mat)
		}
		return(out)
	}
	
	p00max <- omega(x=1,p0_v1=p0_v1,p0_v2=p0_v2)
	p00min <- omega(x=-1,p0_v1=p0_v1,p0_v2=p0_v2)
	
	
	if ((abs(p00-p00max)<tolerance) & (p00>p00max)) p00 <- p00max
	if ((abs(p00-p00min)<tolerance) & (p00<p00min)) p00 <- p00min
	
	if ((p00<p00min) | (p00>p00max) | is.na(p00) ) {
		names(p0_v1) <- "p0_v1"
		names(p0_v2) <- "p0_v2"
        names(p00) <- "p00"
		names(p00max) <- "p00max"
		names(p00max) <- "p00min"
		message(p0_v1) 
		message(p0_v2)
		message(p00)
		message(p00max)
		message(p00min)
		message(.Machine$double.eps)
		stop("Error in omega_inv: p00 out of bounds!!")
		#stop()
	}

	out <- uniroot(omega_root,p0_v1=p0_v1,p0_v2=p0_v2,p00=p00,correlation=correlation,interval=interval,...)
	
	if (only.value) out <- out$root
	return(out)
}