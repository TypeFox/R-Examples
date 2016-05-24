qCompound <-
function(p,parent,compound,compoundDist,params,...) {
    if (!exists(paste("p",parent,sep=""))) {
        return(paste("The parent distribution",parent,"doesn't exist"))
    }
    if (!is.element(compound,compoundDist)) {
        return(paste("The discrete distribution",compound,"doesn't exist"))
    }

    l <- p[p<0|p>1]
    if (length(l)>0) stop("Parameter p is probability")

    xval <- double(length(p))
    Finv <- get(paste("q", parent, sep = ""), mode = "function")
    phi <- get(paste("pgf", compound, sep = ""), mode = "function")
    phiInv <- get(paste("pgfI", compound, sep = ""), mode = "function")
    xval <- Finv(1-phiInv(1-p*(1-phi(0,params)),params),...)
    return(xval)
}
