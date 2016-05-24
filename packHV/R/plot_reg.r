#' Plot points with the corresponding linear regression line
#'
#' Plots points with the corresponding linear regression line
#'
#' @param x numeric vector
#' @param y numeric vector
#' @param pch type of points
#' @param xlab character string, label of the x axis, \code{NULL} by default
#' @param ylab character string, label of the y axis, \code{NULL} by default
#' @param \dots other arguments to be passed in \code{\link{plot}}
#' @return None
#' @author Hugo Varet
#' @examples
#' plot_reg(cgd$age,cgd$height,xlab="Age (years)",ylab="Height")

# last updated: oct 6, 2012

plot_reg=function(x,y,pch=19,xlab=NULL,ylab=NULL,...){
  # x and y are the two variables (missing values are not supported)
  # ... are arguments to be passed in plot() 
  abs <- deparse(substitute(x))
  ord <- deparse(substitute(y))
  
#  eval(parse(text=paste(abs,"=x",sep="")))
#  eval(parse(text=paste(ord,"=y",sep="")))  
#  eval(parse(text=paste("fit=lm(",ord,"~",abs,")",sep="")))
  
  fit=lm(y~x)
  print(summary(fit))
  if (is.null(xlab)){xlab=abs}
  if (is.null(ylab)){ylab=ord}
  d=data.frame(abs=x,ord=y)
  plot(d$abs,d$ord,pch=pch,xlab=xlab,ylab=ylab,...)
  abline(fit,lwd=2)
}

#x1=runif(20)
#y1=0.5+1*x1+rnorm(20,sd=0.1)
#plot_reg(x=x1,y=y1,main="essai",xlab="axe des abscisses",ylab="axe des ordonnées")
#plot_reg(x=x1,y=y1,main="essai")
#y=y1
#x=x1
#plot_reg(y,x,main="essai")
#data=data.frame(abs=x,ord=y)
#plot_reg(data[,"abs"],data[,"ord"],main="essai")
#
#data=data.frame(x,y)
#plot_reg(data$y,data$x,xlab="data$x",ylab="data$y")
#plot_reg(data$y,data$x)
