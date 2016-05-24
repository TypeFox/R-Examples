
NULL

#'
#'  Converts a random variable \code{x} extracted by a population represented by the sample \code{data} or \code{sample}
#'  to a normally-distributed variable with assigned mean and standard deviation or vice versa in case \code{inverse} is \code{TRUE}
#'    
#' @param x value or vector of values to be converted 
#' @param data a sample of data on which a non-parametric probability distribution is estimated
#' @param cpf cumulative probability distribution. If \code{NULL} (default) is calculated as \code{\link{ecdf}(data)}
#' @param mean mean (expected value) of the normalized random variable. Default is 0.
#' @param sd standard deviation of the normalized random variable. Default is 1.
#' @param inverse  logical value. If \code{TRUE} the function works inversely (the opposite way). Default is \code{FALSE}.
#' @param step vector of values in which step discontinuities of the cumulative probability function occur. Default is \code{NULL}
#' @param prec amplitude of the neighbourhood of the step discontinuities where cumulative probability function is treated as non-continuous.
#' @param type see \code{\link{quantile}}
#' @param extremes logical variable. 
#'  If \code{TRUE} (default) the probability or frequency is multiplied by \deqn{\frac{N}{N+1}} where \eqn{N} is the length of \code{data}
#' @param sample a character string or \code{NULL} containing sample or probability distribution information. 
#' Default is \code{NULL} 
#' 
#' @export 
#' 
#' 
#' @author Emanuele Cordano, Emanuele Eccel
#' @return the normalized variable or its inverse   
#'   
#'     
#'  @note This function makes a Marginal Gaussianization. See the R code for further details


#converts a random variable distributed according to data to normally distributed variable





normalizeGaussian <-
function(x=0,data=x,cpf=NULL,mean=0,sd=1,inverse=FALSE,step=NULL,prec=10^-4,type=3,extremes=TRUE,sample=NULL) {
	
	
	
	if (extremes) {
		f=length(data)/(length(data)+1)
	} else {
		f=1
	}
	
	
	
	if (inverse) {
	
		out <- x*NA
		
		qx <- pnorm(x[!is.na(x)])/f # check the extermes!! extremes are 
		
		qx[qx>1] <- 1
		
		out[!is.na(x)] <- quantile(x=data,probs=qx,na.rm=TRUE,names=FALSE,type=type)
		
		
		
		
	} else {
		
		
		
		
		
		if (is.null(cpf)) cpf <- ecdf(data)
		qx <- cpf(x)*f
		# > spline(x=x,y=e(x),xout=y)$y
		
		if (!is.null(step)) {
			for (s in step) {
				for (i in 1:length(x)) {
					
					if (!is.na(x[i])) {
						if (abs(x[i]-s)<prec){
							
							qx2 <-cpf(s+prec)
							qx1 <- cpf(s-prec)
							qx[i] <- runif(1,min=qx1,max=qx2)
							
						}
					}
				}
			}
			
		}
		
		
		
		out <- qnorm(qx,mean=mean,sd=sd)
		
		
		
	}

#	names(out) <- names(x)
	
	
	return(out)
	
}

