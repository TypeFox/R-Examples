#'print.EstTime
#'
#'Generates output for EstTime object, see: \code{\link{EstTime}}
#'
#'
#'@param x a EstTime object, printing it will result in a frequency plot. See help(EstTime)
#'@param pos position of legend, same arguments avaible as in the argument pos from the plot() function.
#'@param ... further arguments passed to or from other methods.
#'@export


print.EstTime<-function(x, pos=NA, ...){
pos<-attributes(x)$pos

plot(smooth(x[[2]],"3RS3R")~x[[1]], typ="l", ylim=c(0,1), xlim=c(min(x[[1]]), max(x[[1]])), xlab="number of time points", ylab="expected")

lines(smooth(x[[3]],"3RS3R")~x[[1]],lty=2)
legend(pos, "cell problems", c("zero frequencies", "low frequencies"), lty=c(1,2))
title("Expected number of cases with cell frequency problems")
}



