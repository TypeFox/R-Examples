gnmfLevel3Parse <- function(HW,margin,nRows,nCols,k,numChains) {
    # If there is only one repeat, just return the input,
    # no parsing necessary.
    if ( numChains == 1 ) {
        return(HW)
    }

    # Parse along rows (margin==1) or columns (margin==2)
    for ( i in 1:numChains ) {
        # Obtain submatrix from the input matrix.
        if ( margin == 1 ) {
            subMat <- HW[((i-1)*k+1):(i*k),]
        } else if ( margin == 2 ) {
            subMat <- HW[,((i-1)*k+1):(i*k)]
        } else {
            print("gnmfLevel3Parse: margin must be either 1 or 2!")
            return(NA)
        }

        # Append to output list.
        if ( i == 1 ) {
            returnVal <- subMat
#            nameList  <- sprintf("%d",i)
        } else if ( i == 2 ) {
            returnVal <- list(returnVal,subMat)
#            nameList  <- list(nameList,sprintf("%d",i))
        } else {
            returnVal[[i]] <- subMat
#            nameList[[i]]  <- sprintf("%d",i)
        }
    }
#    names(returnVal) <- nameList
    return(returnVal)
}
