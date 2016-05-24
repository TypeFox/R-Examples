ikcirt.rndData1 <-
function(N, qTypes, mxDelta, ns) {
    
    mxData <- matrix(NA, sum(ns), N)
    for(i in 1:N) {
        aa <- NULL
        for(j in 1:length(ns)) {
            if(qTypes[j] == "R") {
                this.rank <- sample(I(1:ns[j]))
            }
            if(qTypes[j] == "M") {
                this.rank <- sample( c(1, ns[j], rep(NA, ns[j]-2)) )
            }
            
            aa <- c( aa, this.rank )
        }
        mxData[ , i] <- aa
    }
    mxData
    
    
    return(mxData)
}
