#' Plot for an epplabOutlier Object
#' 
#' Visualizes which observations are considered as outliers and how often for
#' an \code{epplabOutlier} object.
#' 
#' 
#' @name plot.epplabOutlier
#' @aliases plot.epplabOutlier plot,epplabOutlier-method
#' @docType methods
#' @param x Object of class \code{epplabOutlier}.
#' @param col Which colors should be used for non-outliers and outliers.
#' Default is white and black.
#' @param outlier Logical if only observations considered as outliers at least
#' once should be plotted or all. Default is \code{TRUE}.
#' xlab Possible x axis label. Default is none.
#' ylab Possible x axis label. Default is none.
#' @param xlab Sets the x-axis label.
#' @param ylab Sets the y-axis label.
#' @param ... Graphical parameters passed on to \code{\link{image}}.
#' @author Daniel Fischer and Klaus Nordhausen
#' @seealso \code{\link{EPPlabOutlier}}, \code{\link{image}}
#' @keywords methods hplot
#' @examples
#' 
#' # creating data with 3 outliers
#' n <-300 
#' p <- 10
#' X <- matrix(rnorm(n*p),ncol=p)
#' X[1,1] <- 9
#' X[2,4] <- 7 
#' X[3,6] <- 8
#' # giving the data rownames, obs.1, obs.2 and obs.3 are the outliers.
#' rownames(X) <- paste("obs",1:n,sep=".")
#' 
#' PP<-EPPlab(X,PPalg="PSO",PPindex="KurtosisMax",n.simu=20, maxiter=20)
#' OUT<-EPPlabOutlier(PP, k = 3, location = median, scale = mad)
#' plot(OUT)
#' 
#' @export
plot.epplabOutlier <- function(x,col=c("white","black"),outlier=TRUE, xlab="", ylab="", ...) 
{
  # If only outlier are requested, then the outlier matrix is limitied to those, where
  # the rowsum is larger than 0.
  ifelse(outlier, X <- x$outlier[apply(x$outlier,1,sum)>0,,drop=FALSE] ,  X <- x$outlier)
 
  p <- ncol(X)
  n <- nrow(X)
  if (n==0) stop('Option "outlier=TRUE" not meaningful since no outliers were detected')
  X2 <- t(X[n:1,, drop=FALSE])

  image(1:p,1:n,X2, col=col, axes=FALSE, xlab=xlab, ylab=ylab,...)
  box()
  axis(2, labels=colnames(X2), at=1:n,...)
  axis(1, labels=rownames(X2), at=1:p,...)    
}
