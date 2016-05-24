#' mypar
#' 
#' Called without arguments, this function optimizes graphical parameters
#' for the RStudio plot window. \code{bigpar} uses big fonts which are good for presentations.
#' 
#' @aliases bigpar
#' 
#' @param a the first entry of the vector passed to \code{mar}
#' @param b the second entry of the vector passed to \code{mar}
#' @param brewer.n parameter \code{n} passed to \code{\link{brewer.pal}}
#' @param brewer.name parameters \code{name} passed to \code{\link{brewer.pal}}
#' @param cex.lab passed on to \code{\link{par}}
#' @param cex.main passed on to \code{\link{par}}
#' @param cex.axis passed on to \code{\link{par}}
#' @param mar passed on to \code{\link{par}}
#' @param mgp passed on to \code{\link{par}}
#' @param ... other parameters passed on to \code{\link{par}}
#' 
#' @author Rafael A. Irizarry
#' 
#' @examples
#' mypar()
#' plot(cars)
#' bigpar()
#' plot(cars)

mypar <- function(a=1,b=1,brewer.n=8,brewer.name="Dark2",cex.lab=1,cex.main=1.2,cex.axis=1,mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0),...){
    par(mar=mar,mgp=mgp,cex.lab=cex.lab,cex.main=cex.main,cex.axis=cex.axis)
    par(mfrow=c(a,b),...)
    palette(RColorBrewer::brewer.pal(brewer.n,brewer.name))
}

bigpar <- function(a=1,b=1,brewer.n=8,brewer.name="Dark2",cex.lab=2,cex.main=2,cex.axis=1.5,mar=c(5.1,5.1,3.5,2.1),mgp=c(3,1,0),...){
    par(mar=mar,mgp=mgp,cex.lab=cex.lab,cex.main=cex.main,cex.axis=cex.axis)
    par(mfrow=c(a,b),...)
    palette(RColorBrewer::brewer.pal(brewer.n,brewer.name))
}
 
