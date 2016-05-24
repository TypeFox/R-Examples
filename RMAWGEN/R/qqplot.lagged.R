NULL
#' 
#' This function creates a Q-Q plot of  the  \code{lag}-lag moving cumulative addition of the values in the samples \code{x,y,z} 
#' 
#' @param x,y samples. If \code{x} is a data frame, \code{y} and \code{z} can be omitted.
#' @param z further samples organized as a list
#' @param when (integer) inidices of \code{x} and \code{y} on which the Q-Q plot is made. 
#' @param lag lag (current index included) on whose value  the addition is made.
#' @param pch a vector of plotting characters or symbols: see \code{\link{points}}
#' @param ... further arguments for \code{\link{qqplot}}
#' 
#' @return the Q-Q plot
#' @export
#' @seealso \code{\link{qqplot}}

# THIS FUCTION IS TO CHECK??? 

qqplot.lagged <- function(x=rnorm(1000),y=rnorm(1000),z=NULL,when=1:length(x),lag=1,pch=1,...){
	
	if (lag<=0) {
		
		stop("Error in qqplot.lagged: lag less than 1!!")
		
	}
	
	if (is.data.frame(x)) {
		
		xdata <- x
		x <- xdata[,1]
		y <- xdata[,2]
		
		if (ncol(xdata)>2) {
			
			z <- list()
			for (i in 1:(ncol(xdata)-2)) {
				
				z[[i]] <- xdata[,i+2]
			}
		}
		
		
	}
	
	
	
	
	xl <- array(0,length(x)-lag+1)
	yl <- array(0,length(y)-lag+1)
	
	if (is.list(z)) {
		pch <- 1:(length(z)+1)
		zl <- list() 
		for (i in 1:length(z)) {
			
			zl[[i]] <- array(0,length(z[[i]])-lag+1)
		}
		
	}
	for (l in 0:(lag-1)) {
		
		xl <- xl+x[(lag-l):(length(x)-l)]
		yl <- yl+y[(lag-l):(length(y)-l)]
		
		if (is.list(z)) {
			
			for (i in 1:length(z)) {
				
				
				zl[[i]] <- zl[[i]]+z[[i]][(lag-l):(length(z[[i]])-l)]
				
			}
		}
		
	} 
	
	whenl <- when[when>=lag]-lag+1
	out <- qqplot(xl[whenl],yl[whenl],pch=pch[1],...)
	
	if (is.list(z)) {
		
		xs <- sort(xl[whenl])
		for (i in 1:length(z)) {
			
			zs <- sort(zl[[i]][whenl])
			points(xs,zs,pch=pch[i+1],...)
			
		}
		
		
	}
	
	return(0)
	
	
}




#qqplot(generation00$prec_mes[iseason,istation],generation00$prec_gen[iseason,istation],main=title,xlab="observed",ylab="generated",cex=cex,pch=pch)