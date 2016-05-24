pCompound <-
function(q,parent,compound,compoundDist,params,...) {
    if (!exists(paste("p",parent,sep=""))) {
        return(paste("The parent distribution",parent,"doesn't exist"))
    }
    if (!is.element(compound,compoundDist)) {
        return(paste("The discrete distribution",compound,"doesn't exist"))
    }
    xval <- double(length(q))
    F <- get(paste("p",parent,sep=""), mode = "function")
    phi <- get(paste("pgf",compound,sep=""), mode = "function")
    xval <- (1-phi(1-F(q,...),params))/(1-phi(0,params))
    return(xval)
}
