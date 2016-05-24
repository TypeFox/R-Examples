# Function name: acceptanceProbability
# Author: Yang Li <yang.li@rug.nl>
# Maintainer: Yang Li <yang.li@rug.nl>
# Version: 1.0.0
# Date: 10 Dec. 2007


acceptanceProbability <- function( designScore, newDesignScore, method,
                                   temperature )
{   
    #used by designGG function (Simulated annealing process)
    if( method == "MH" )
    {
        acc.prob <- newDesignScore/designScore
    }
    else
    {
        if( method == "SA" )
        {
            acc.prob <- (newDesignScore/designScore)^(1/temperature)
        }
        else
        {
            stop("This acceptance probability method has not been implemented")
        }
    }
    
    if( is.na(acc.prob) )
      acc.prob <- -1

    return( acc.prob )
}