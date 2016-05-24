hCompound <-
function(x, parent, compound,compoundDist,params, ...) {
    if (!exists(paste("p",parent,sep=""))) {
        return(paste("The parent distribution",parent,"doesn't exist"))
    }
    if (!is.element(compound,compoundDist)) {
        return(paste("The discrete distribution",compound,"doesn't exist"))
    }
    dCompound(x,parent,compound,compoundDist,params,...)/(1-pCompound(x,parent,compound,compoundDist,params,...))
}
