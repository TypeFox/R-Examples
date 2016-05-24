# Function name: pairLevel
# Author: Yang Li <yang.li@rug.nl>
# Maintainer: Yang Li <yang.li@rug.nl>
# Version: 1.0.0
# Date: 10 Dec. 2007


pairLevel <- function( xxx, rilNames )
{
    #This is used by designScore function
    xxx         <- as.matrix(xxx)
    pair.level  <- NULL
    rowNames    <- NULL
    for( i in 1:(nrow(xxx)/2) )
    {
        pair.level <- rbind( pair.level, xxx[(2*i-1),] - xxx[2*i,] )
        rowNames   <- c( rowNames, paste(rilNames[(2*i-1)], "&",
                         rilNames[(2*i)], sep="") )
    }
    rownames(pair.level) <- rowNames
    return( pair.level )
}