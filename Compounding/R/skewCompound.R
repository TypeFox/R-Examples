skewCompound <-
function(parent,compound,compoundDist,params,...) {
    if (!exists(paste("p",parent,sep=""))) {
        return(paste("The parent distribution",parent,"doesn't exist"))
    }
    if (!is.element(compound,compoundDist)) {
        return(paste("The discrete distribution",compound,"doesn't exist"))
    }
   m1 <- meanCompound(parent,compound,compoundDist,params,...)
   m3 <- momentCompound(3,parent,compound,compoundDist,params,...)
   sig2 <- varCompound(parent,compound,compoundDist,params,...)
   return((m3-3*m1*sig2-m1^3)/sig2^(3/2))
}
