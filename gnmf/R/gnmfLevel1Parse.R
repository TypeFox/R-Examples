gnmfLevel1Parse <- function(HW,margin,nRows,nCols,alphas,ranks,numChains) {
    # Parse along rows (margin==1) or columns (margin==2)
    numAlphas      <- length(alphas)
    if ( length(ranks) == 1 ) {
        alphaBlockSize <- ranks*numChains
    } else {
        alphaBlockSize <- sum(ranks[1]:ranks[2])*numChains
    }

    for ( i in 1:numAlphas) {
        # Determine offset into rows/columns due to current rank.
        a <- alphas[i]
        if ( i == 1 ) {
            alphaOffset <- 0
        } else {
            alphaOffset <- (i-1)*alphaBlockSize
        }

        # Obtain submatrix from the input matrix.
        if ( margin == 1 ) {
            subMat <- HW[(alphaOffset+1):(alphaOffset+alphaBlockSize),]
        } else if ( margin == 2 ) {
            subMat <- HW[,(alphaOffset+1):(alphaOffset+alphaBlockSize)]
        } else {
            print("gnmfLevel1Parse: margin must be either 1 or 2!")
            return(NA)
        }
        parsedMatList <- gnmfLevel2Parse(subMat,margin,nRows,nCols,ranks,numChains)

        # Append to output list and maintain list of names.
        if ( i == 1 ) {
            returnVal <- parsedMatList

            # If there is only one rank, just return the input,
            # no further parsing necessary.
            if ( numAlphas == 1 ) {
                return(returnVal)
            }

            # Otherwise maintain a name list for parsing.
            nameList  <- sprintf("%g",a)
        } else if ( i == 2 ) {
            returnVal <- list(returnVal,parsedMatList)
            nameList  <- list(nameList,sprintf("%g",a))
        } else {
            returnVal[[i]] <- parsedMatList
            nameList[[i]]  <- sprintf("%g",a)
        }
    }
    names(returnVal) <- nameList
    return(returnVal)
}
