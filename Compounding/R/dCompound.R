dCompound <-
function(x,parent,compound,compoundDist,params,...) {
    if (!exists(paste("p",parent,sep=""))) {
        return(paste("The parent distribution",parent,"doesn't exist"))
    }
    if (!is.element(compound,compoundDist)) {
        return(paste("The discrete distribution",compound,"doesn't exist"))
    }
    xval <- double(length(x))
    f <- get(paste("d", parent, sep = ""), mode = "function")
    F <- get(paste("p", parent, sep = ""), mode = "function")
    phi <- get(paste("pgf",compound,sep=""), mode = "function")
    phiD <- get(paste("pgfD",compound,sep=""), mode = "function")
    xval <- phiD(1-F(x,...),params)*f(x,...)/(1-phi(0,params))
    return(xval)
}
