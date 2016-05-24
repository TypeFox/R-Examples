rCompound <-
function(n, parent, compound,compoundDist,params, ...) {
    if (!exists(paste("p",parent,sep=""))) {
        return(paste("The parent distribution",parent,"doesn't exist"))
    }
    if (!is.element(compound,compoundDist)) {
        return(paste("The discrete distribution",compound,"doesn't exist"))
    }
    
    if (n<0) stop("Parameter n must be positive")
    if(!(abs(n-round(n)) < .Machine$double.eps^0.5)) stop("Parameter n must be positive integer")
    zval <- runif(n)
    xval <- qCompound(zval,parent,compound,compoundDist,params,...)
    return(xval)
}
