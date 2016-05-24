#' Plot chosen coefficients of frequenct domain operator. Draws all pairs \eqn{ \{ (i,j) | i \in x, j\in y }.
#'
#' @title Plot a frequency domain operator
#' @param x Frequency domain object
#' @param xcoef vector of horizontal coefficients to plot
#' @param ycoef vector of vertical coefficients to plot
#' @param type 'cartesian' or 'polar' coordinates
#' @param ... other parameters to \code{\link{plot}}
#' @importFrom graphics lines par plot
#' @export
plot.freqdom = function(x, xcoef=NULL, ycoef=NULL, type='cartesian', ...){
  SD = x

  if (!is.freqdom(SD))
    stop("SD must be a freqdom object")

  D = dim(SD$operators)
  if (is.null(xcoef))
    xcoef = 1:D[2]
  if (is.null(ycoef))
    ycoef = 1:D[3]

  if (!is.vector(xcoef) || !is.numeric(xcoef))
    stop("x must be a numeric vector")
  if (!is.vector(ycoef) || !is.numeric(ycoef))
    stop("x must be a numeric vector")
  
  ylim = min(Re(SD$operators),Im(SD$operators))
  ylim = c(ylim,max(Re(SD$operators),Im(SD$operators)))
  ylim = max(abs(ylim))
  ylim = c(-ylim,ylim)

  par(mfrow=c(length(xcoef),length(ycoef)))
  for (i in xcoef){
    for (j in ycoef){
      plot.fd.coef(SD,i,j,ylim=ylim,type=type)
    }
  }
}

plot.fd.coef = function(SD, x, y, ylim=NULL, type='cartesian', ...){
  if (type=='cartesian'){
    f1 = Re
    f2 = Im
  }
  if (type=='polar'){
    f1 = Mod
    f2 = Arg
  }
  if (is.null(ylim)){
    ylim = min(f1(SD$operators[,x,y]),f2(SD$operators[,x,y]))
    ylim = c(ylim,max(f1(SD$operators[,x,y]),f2(SD$operators[,x,y])))
  }
  plot(SD$freq, f1(SD$operators[,x,y]),t='l',xlab="frequency",ylab="value",ylim=ylim, ...)
  lines(SD$freq, f2(SD$operators[,x,y]),t='l',col=2)
  title(paste("Coeficient (",x,",",y,") of Frequency Domain Operator", sep=""))
}
