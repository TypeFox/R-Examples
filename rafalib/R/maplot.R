#' Bland Altman plot aka MA plot
#'
#' Takes two vectors x and y and plots 
#' M=y-x versus A=(x+y)/2. 
#' If the vectors a more longer than length n the data is sampled to size n. 
#' A smooth curve is added to show trends.
#' 
#' @param x a numeric vector
#' @param y a numeric vector
#' @param n a numeric value. If \code{length(x)} is larger than \code{n}, the \code{x} and \code{y} are sampled down.
#' @param subset index of the points to be plotted 
#' @param xlab a title for the x axis
#' @param ylab a title for the y axis
#' @param curve.add if \code{TRUE} a smooth curve is fit to the data and displayed. The function \code{\link{loess}} is used to fit the curve.
#' @param curve.col a numeric value that determines the color of the smooth curve
#' @param curve.span is passed on to \code{\link{loess}} as the \code{span} argument
#' @param curve.lwd the line width for the smooth curve
#' @param curve.n a numeric value that determines the sample size used to fit the curve. This makes fitting the curve faster with large datasets
#' @param ... further arguments passed to \code{\link{plot}}
#' 
#' @author Rafael A. Irizarry
#'  
#' @examples
#' n <- 10000
#' signal <- runif(n,4,15)
#' bias <- (signal/5 - 2)^2 
#' x <- signal + rnorm(n)
#' y <- signal + bias + rnorm(n)  
#' maplot(x,y)

maplot <- function(x,y,n=10000,subset=NULL,xlab=NULL,ylab=NULL,
                   curve.add=TRUE,curve.col=2,curve.span=1/2,
                   curve.lwd=2,curve.n=2000,...){
    if(length(x)!=length(y)) stop("Length of 'x' and 'y' must be the same.")  
    if(is.null(xlab)) xlab="A"
    if(is.null(ylab)) ylab="M"
    if(!is.null(subset)){
        x=x[subset]
        y=y[subset]
    }
    if(length(x)>n){
        ind=sample(length(x),10000)
        x=x[ind]
        y=y[ind]
    }
    m=y-x
    a=(x+y)/2
    plot(a,m,xlab=xlab,ylab=ylab,...)
    if(curve.add){
        o=order(a)
        aa=a[o]
        mm=m[o]
        o=seq(1,length(aa),len=curve.n)
        aa=aa[o]
        mm=mm[o]
        fit1=loess(mm~aa,span=curve.span,degree=1)
        lines(aa,fit1$fitted,col=curve.col,lwd=curve.lwd)
    }
}
