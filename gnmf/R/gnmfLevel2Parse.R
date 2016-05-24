gnmfLevel2Parse <- function(HW,margin,nRows,nCols,ranks,numChains) {
    # Parse along rows (margin==1) or columns (margin==2)
    if ( length(ranks) == 1 ) {
        rankRange <- ranks
    } else {
        rankRange <- ( ranks[1]:ranks[2] )
    }
    numRanks  <- length(rankRange)
    for ( i in 1:numRanks ) {
        # Determine offset into rows/columns due to current rank.
        r <- rankRange[i]
        if ( i == 1 ) {
            rankOffset <- 0
        } else {
            rankOffset <- (r-1)*numChains
        }

        # Obtain submatrix from the input matrix.
        if ( margin == 1 ) {
            subMat <- HW[(rankOffset+1):(rankOffset+(r*numChains)),]
        } else if ( margin == 2 ) {
            subMat <- HW[,(rankOffset+1):(rankOffset+(r*numChains))]
        } else {
            print("gnmfLevel3Parse: margin must be either 1 or 2!")
            return(NA)
        }
        parsedMatList <- gnmfLevel3Parse(subMat,margin,nRows,nCols,r,numChains)

        # Append to output list and maintain list of names.
        if ( i == 1 ) {
            returnVal <- parsedMatList

            # If there is only one rank, just return the input,
            # no further parsing necessary.
            if ( numRanks == 1 ) {
                return(returnVal)
            }

            # Otherwise maintain a name list for parsing.
            nameList  <- sprintf("%d",r)
        } else if ( i == 2 ) {
            returnVal <- list(returnVal,parsedMatList)
            nameList  <- list(nameList,sprintf("%d",r))
        } else {
            returnVal[[i]] <- parsedMatList
            nameList[[i]]  <- sprintf("%d",r)
        }
    }
    names(returnVal) <- nameList
    return(returnVal)
}
