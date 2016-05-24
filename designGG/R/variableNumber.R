# Function name: variableNumber 
# Author: Yang Li <yang.li@rug.nl>
# Maintainer: Yang Li <yang.li@rug.nl>
# Version: 1.0.0
# Date: 10 Dec. 2007


variableNumber <- function( nEnvFactors )
{
    #This is used by designGG function
    
    n.var <- 0

    for( i in 1:(nEnvFactors+1) )
    {
      n.var <- n.var + choose(nEnvFactors+1,i)
    }
    return( n.var )
}