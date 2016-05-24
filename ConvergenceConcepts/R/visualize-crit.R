 # To draw the criterion curve
visualize.crit <- function(crit,plotfunc=plot,...) {

  plotfunc(crit,type="l",xlab="n",ylab='',lwd=1,...)

}
