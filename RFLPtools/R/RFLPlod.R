###############################################################################
## Remove bands below LOD
###############################################################################

RFLPlod <- function(x, LOD){
    if(length(LOD) > 1){
        warning("Only first element of 'LOD' is used.")
        LOD <- LOD[1]
    }
    
    message(paste(sum(x$MW < LOD), "bands were removed."))
    x[x$MW >= LOD, ]
}
