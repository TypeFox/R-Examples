plot.PVE <- function(x,xlab="Number of Features",ylab="PVE",...) {

    if (!inherits(x,"PVE")) {
        stop("'x' must be of class 'PVE'")
    }

    plot(x=x$J,y=x$PVEs,xlab=xlab,ylab=ylab,...)

}
