a2sf <- function(x,tt,xlim,npts=1000) {
#
# Convert a function built by approx() to a step function.
#
if(missing(tt)) {
	if(missing(xlim))
		stop("One of tt and xlim must be specified.\n")
	tt <- seq(xlim[1],xlim[2],length=npts)
}
y <- x(tt)
dy <- diff(y)
iii <- which(dy!=0)
if(length(iii)==0) {
	knots <- range(tt)
	values <- c(y[1],y[1],y[1])
} else {
	knots <- tt[iii]
	values <- c(y[1],y[1+iii])
}
stepfun(x=knots,y=values)
}
