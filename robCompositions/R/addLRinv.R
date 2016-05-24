#' Additive logistic transformaton
#' 
#' Inverse additive log-ratio transformation, often called additive logistic
#' transformation.
#' 
#' The function allows also to preserve absolute values when class info is
#' provided. Otherwise only the relative information is preserved.
#' 
#' @param x data set, object of class \dQuote{alr}, \dQuote{matrix} or
#' \dQuote{data.frame}
#' @param cnames column names. If the object is of class \dQuote{alr} the
#' column names are chosen from therein.
#' @param ivar index of the rationing part. If the object is of class
#' \dQuote{alr} the column names are chosen from therein. If not and ivar is
#' not provided by the user, it is assumed that the rationing part was the last
#' column of the data in the simplex.
#' @param useClassInfo if FALSE, the class information of object \code{x} is
#' not used.
#' @return the transformed data matrix
#' @export
#' @author Matthias Templ
#' @seealso \code{\link{isomLRinv}}, \code{\link{cenLRinv}},
#' \code{\link{cenLR}}, \code{\link{addLR}}
#' @references Aitchison, J. (1986) \emph{The Statistical Analysis of
#' Compositional Data} Monographs on Statistics and Applied Probability.
#' Chapman \& Hall Ltd., London (UK). 416p.
#' @keywords manip
#' @examples
#' 
#' data(arcticLake)
#' x <- arcticLake
#' x.alr <- addLR(x, 2)
#' y <- addLRinv(x.alr)
#' ## This exactly fulfills:
#' addLRinv(addLR(x, 3))
#' data(expenditures)
#' x <- expenditures
#' y <- addLRinv(addLR(x, 5))
#' head(x)
#' head(y)
#' ## --> absolute values are preserved as well.
#' 
#' ## preserve only the ratios:
#' addLRinv(x.alr, ivar=2, useClassInfo=FALSE)
#' 
#' 
#' 
addLRinv <- function(x, cnames=NULL, ivar=NULL, useClassInfo=TRUE){
	if(class(x) == "alr" & useClassInfo==TRUE){
		
		xalr <- x$x.alr
		ivar <- x$ivar
		dat <- exp(xalr)*x$varx
		## correct order of the variables:
		if(ivar == dim(xalr)[2]+1){ 
			dat <- cbind(dat, x$varx)
		} else if(ivar == 1){  
			dat <- cbind(x$varx, dat)			
		} else{
			dat <- cbind(dat[,1:(ivar-1)], x$varx, dat[,(ivar):(dim(xalr)[2])])
		}
		colnames(dat) <- x$cnames
		if(class(x$x.alr) == "data.frame") dat <- as.data.frame(dat)
	} else if(class(x)=="alr" & useClassInfo == FALSE){
		if(is.null(ivar)) stop("object ivar must be provided \n because object x is not from class alr")
		xalr <- x$x.alr
		#if(is.null(cnames)) cnames <- c(colnames(x), "rat")
		#if(length(cnames)==1) cnames <- paste("V", 1:dim(x)[2]+1, sep="")	
		#if(length(cnames) != dim(x)[2] + 1 | length(cnames) < 2) stop(paste("cnames must be of length", dim(x)[2]+1))
		rat <- rowSums(exp(xalr)) + 1 
		dat <- exp(xalr)/rat
		## correct order of the variables:
		if(ivar == dim(xalr)[2]+1){ 
			dat <- cbind(dat, 1/rat)
		} else if(ivar == 1){
			dat <- cbind(1/rat, dat)			
		} else{
			dat <- cbind(dat[,1:(ivar-1)], 1/rat, dat[,(ivar):(dim(xalr)[2])])
		}
		#colnames(dat) <- x$cnames
	} else if(class(x) != "alr"){
		if(dim(x)[2] < 2) stop("data must be of dimension greater equal 2")
		#if(useClassInfo) warning("x is not from class alr, absolute values are not preserved and column names may not be respected")
		if(is.null(ivar)){ 
			warning(paste("object ivar is not provided \n it is assigned to ", dim(x)[2]+1, sep=""))
			ivar <- dim(x)[2]+1
		}
		xalr <- x
		if(is.null(cnames)) cnames <- c(colnames(x), "rat")
		if(length(cnames)==1) cnames <- paste("V", 1:dim(x)[2]+1, sep="")	
		if(length(cnames) != dim(x)[2] + 1 | length(cnames) < 2) stop(paste("cnames must be of length", dim(x)[2]+1))
		rat <- rowSums(exp(xalr)) + 1 
		dat <- exp(xalr)/rat
		## correct order of the variables:
		if(ivar == dim(xalr)[2]+1){ 
			dat <- cbind(dat, 1/rat)
		} else if(ivar == 1){
			dat <- cbind(1/rat, dat)			
		} else{
			dat <- cbind(dat[,1:(ivar-1)], 1/rat, dat[,(ivar):(dim(xalr)[2])])
		}		
		colnames(dat) <- cnames
	}
	
	
	
	return(dat)
}

