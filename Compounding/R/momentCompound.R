momentCompound <-
function(k, parent, compound,compoundDist, params, ...) {
    if (!exists(paste("p",parent,sep=""))) {
        return(paste("The parent distribution",parent,"doesn't exist"))
    }
    if (!is.element(compound,compoundDist)) {
        return(paste("The discrete distribution",compound,"doesn't exist"))
    }
    if (k<0) stop("Parameter k must be positive")
    if(!(abs(k-round(k))<.Machine$double.eps^0.5)) stop("Parameter k must be positive integer")

    Finv <- get(paste("q", parent, sep = ""), mode = "function")
    phi <- get(paste("pgf",compound,sep=""), mode = "function")
    phiD <- get(paste("pgfD",compound,sep=""), mode = "function")
    fint <- function(x) phiD(1-x,params)*(Finv(x,...))^k/(1-phi(0,params))
    return(integrate(fint,lower=0,upper=1)$value)
}
