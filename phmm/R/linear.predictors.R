#' PHMM Design
#' 
#' Internal function for extracting the linear predictors of the PHMM model
#' from an object of class \code{phmm} returned by \code{phmm}.
#' 
#' 
#' @param x an object of class \code{phmm}.
#' @return A vector of estimates from call to \code{phmm} of
#' \eqn{\beta'x_{ij}+w_{ij}'b_{i}}{beta'x[ij]+w[ij]'b[i]}.
#' @seealso \code{\link{phmm}}
#' @keywords survival
linear.predictors <- function (x) UseMethod("linear.predictors")
linear.predictors.phmm <- function(x){
    #Function to compute linear predictors
    #cluster is assumed to be sorted
	
	z=x$Z; beta=x$coef; w=x$W; b=as.matrix(x$bhat.long)
	
    z <- as.matrix(z)
    wb <- matrix(0,nrow=nrow(w),ncol=ncol(w))
    for(i in 1:ncol(w)){
	  if(length(w[,i])!=length(b[,i])) stop("length(w[,i])!=length(b[,i])")
	  wb[,i] <- w[,i]*b[,i]
	  }
    if(!is.null(dim(wb))) wb <- apply(wb,1,sum)
    return(z%*%beta+wb)
    }