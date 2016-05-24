varCompound <-
function(parent,compound,compoundDist,params,...) {
    if (!exists(paste("p",parent,sep=""))) {
        return(paste("The parent distribution",parent,"doesn't exist"))
    }
    if (!is.element(compound,compoundDist)) {
        return(paste("The discrete distribution",compound,"doesn't exist"))
    }
    m1 <- momentCompound(1,parent,compound,compoundDist,params,...)
    m2 <- momentCompound(2,parent,compound,compoundDist,params,...)
    return(m2-m1^2)
}
