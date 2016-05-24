ikcirt.Ustar2data <-
function(Ustar, qTypes, mxDelta, ns) {
    
    N <- ncol(Ustar)
    mxData <- matrix(NA, sum(ns), N)
    ndxVec <- c(0, cumsum(ns)) + 1
    for(i in 1:N) {
        aa <- NULL
        for(j in 1:length(ns)) {
            theseVals <- Ustar[ ndxVec[j]:(ndxVec[j+1]-1), i ]
            this.rank <- rank(theseVals)
            
            if(qTypes[j] == "R") {
                
            }
            if(qTypes[j] == "M") {
                this.rank[ which( this.rank > 1 & this.rank < max(this.rank) ) ] <- NA
            }
            ### print(this.rank)
            aa <- c( aa, this.rank )
        }
        mxData[ , i] <- aa
    }
    mxData
    return(mxData)
}
