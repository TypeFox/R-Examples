
NULL

#'
#' 
#' Converts precipitation values to "Gaussinized" normally-distributed values taking into account the probability of no precipitation occurences. values  
#' or vice versa in case \code{inverse} is \code{TRUE}
#'    
#' @param x value or vector of values to be converted 
#' @param data a sample of data on which a non-parametric probability distribution is estimated
#' @param cpf cumulative probability distribution. If \code{NULL} (default) is calculated as \code{\link{ecdf}(data)}
#' @param mean mean (expected value) of the normalized random variable. Default is 0.
#' @param sd standard deviation of the normalized random variable. Default is 1.
#' @param inverse  logical value. If \code{TRUE} the function works inversely (the opposite way). Default is \code{FALSE}.
#' @param qnull probability of no precipitation occurence
#' @param valmin minimum value of precipitation to consider a wet day
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
#' @seealso \code{\link{normalizeGaussian}}     
#' @note In the version 1.2.5 of \pkg{RMAWGEN} This function is deprecated and not used.
# This function makes a Marginal Gaussianization accounting for the days with zero precipitation. See the R code for further details
#' 
#' 
#' 
#' 
#' 
#' @examples 
#' library(RMAWGEN)
#' NDATA <- 1000
#' occurence <- as.logical(runif(NDATA)>0.5)
#' prec <- rexp(NDATA,rate=1/3)
#' prec[!occurence] <- 0
#' valmin <- 0.5 #0.01
#' x <- normalizeGaussian_prec(x=prec,valmin=valmin)
#' prec2 <- normalizeGaussian_prec(x=x,data=prec,valmin=valmin,inverse=TRUE)
#' qqplot(prec,prec2)
#'
#' occurence3 <- as.logical(runif(NDATA)>0.5)
#' prec3 <- rexp(NDATA,rate=1/3)
#' prec3[!occurence3] <- 0
#' x3 <- normalizeGaussian_prec(x=prec3,valmin=valmin) 
#' 
#' qqplot(x,x3)
#' abline(0,1)
#' 

#' 







normalizeGaussian_prec <-
function(x=0,data=x,cpf=NULL,mean=0,sd=1,inverse=FALSE,type=3,extremes=TRUE,sample=NULL,qnull=0,valmin=1.0) {
	
	
	
	
	
	if (is.null(valmin)) {
		qnull=0
		valmin=min(data,na.rm=TRUE)-0.001
	}
	
	
	if (extremes) {
		f=length(data[!is.na(data) & data>valmin])/(length(data[!is.na(data) & data>valmin])+1)
	} else {
		f=1
	}
	
	if ((is.null(qnull))) qnull=length(data[!is.na(data) & data<=valmin])/length(data[!is.na(data)])*f
	
	
	
	out <- array(NA,length(x))
	
	if (inverse) {
	
		qx <- x*NA
		
		qx[!is.na(x)] <- pnorm(x[!is.na(x)])/f # check the extermes!! extremes are 
		
		qx[!is.na(x) & qx>1] <- 1
		
		qx0 <- qx[!is.na(x) & qx>qnull] 
		
		qx00 <- (qx0-qnull)/(1-qnull)
		
		
		out[!is.na(x) & qx>qnull] <- quantile(x=data[!is.na(data) & data>valmin],probs=qx00,na.rm=TRUE,names=FALSE,type=type)
		
		out[!is.na(x) & qx<=qnull] <- valmin
		
	
	} else {
		
		
		x0 <- x[!is.na(x) & x>valmin]
		
		
		if (is.null(cpf)) cpf <- ecdf(data[!is.na(data) & data>valmin])
		qx0 <- cpf(x0)*f*(1-qnull)

	
		
		out[!is.na(x) & x>valmin] <- qnorm(qx0+qnull,mean=mean,sd=sd)
		out[!is.na(x) & x<=valmin] <- qnorm(runif(1,min=0,max=qnull),mean=mean,sd=sd)	

		
	}

#	names(out) <- names(x)
	
	
	return(out)
	
}

