#' altman-bland plot
#' @description
#' plots abs(x-y) against (y+x)/2
#' @param x - input intensities
#' @param y - input intensities
#' @param main plotting parameters
#' @param pch - plot character
#' @param log - should the x y axis be log transformed possible values "x" "y" or "xy"
#' @export
#' @examples
#' x <- seq(1:300)/10
#' x <- x + rnorm(length(x),0,0.5)
#' y <- seq(1:300)/10
#' y <- y + rnorm(length(y),0,0.5)
#' altmanbland(y,x)
#'
altmanbland = function(x,y,main="",pch=".",log=""){
  idx<-apply(cbind(x,y),1,function(x){sum(x==0)==0})
  mean  = (x+y)/2
  absdiff = abs( x-y )
  #scatter.smooth(mean[idx],absdiff[idx],span=1/3,col=2,pch=".",cex=2,xlab="(y+x)/2",ylab="abs(x-y)",main=main,log="xy")
  plot(mean,absdiff,log=log,xlab="(y+x)/2",ylab="abs(x-y)",pch="x",cex=0.5,main=main)
  lines(lowess(mean,absdiff),col=2,lwd=2)
}
