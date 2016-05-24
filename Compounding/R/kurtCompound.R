kurtCompound <-
function(parent,compound,compoundDist,params,...) {
    if (!exists(paste("p",parent,sep=""))) {
        return(paste("The parent distribution",parent,"doesn't exist"))
    }
    if (!is.element(compound,compoundDist)) {
        return(paste("The discrete distribution",compound,"doesn't exist"))
    }
   m1 <- meanCompound(parent,compound,compoundDist,params,...)
   m3 <- momentCompound(3,parent,compound,compoundDist,params,...)
   m4 <- momentCompound(4,parent,compound,compoundDist,params,...)
   sig2 <- varCompound(parent,compound,compoundDist,params,...)
   return((m4-4*m1*m3+6*m1^2*sig2+3*m1^4)/sig2^2)
}
