#' smart plot
#'
#' if n > 10,000, make a random subset of 10,000 and plot. You can also specify
#' a specific subset to plot. If length of subset is larger 
#' than \code{n}, a random sample is still used to reduce data size.
#'
#' @param x the x data
#' @param y the y data
#' @param n the number to subset
#' @param subset explicit subset index (optional). 
#' @param xlab title for the x-axis
#' @param ylab title for the y-axis
#' @param ... further parameters passed on to \code{plot}
#' @examples
#'
#' x <- rnorm(1e5)
#' y <- rnorm(1e5)
#' splot(x,y,pch=16,col=rgb(0,0,0,.25))
#' 
splot <- function(x,y,n=10000,subset=NULL,xlab=NULL,ylab=NULL,...){
    if(is.null(xlab)) xlab=deparse(substitute(x))
    if(is.null(ylab)) ylab=deparse(substitute(y))
    if(!is.null(subset)){
        x=x[subset]
        y=y[subset]
    }
    if(length(x)>n){
        ind=sample(length(x),10000)
        x=x[ind]
        y=y[ind]
    }
    plot(x,y,xlab=xlab,ylab=ylab,...)
}
