assignColors <- function(cols,MAT){
    colorN <- length(cols)
    rangeRR <- range(as.vector(MAT[!is.na(as.vector(MAT))]))
    diffRR <- (rangeRR[2]-rangeRR[1])/colorN
    vMAT <- vRET <- as.vector(MAT)
    
    for(i in 0:(colorN-1)) vRET[vMAT<=(rangeRR[2]-i*diffRR)] <- cols[colorN-i]

    return(t(matrix(vRET,nrow(MAT),ncol(MAT))))
}
