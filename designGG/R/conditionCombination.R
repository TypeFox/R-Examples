# Function name: conditionCombination
# Author: Yang Li <yang.li@rug.nl>
# Maintainer: Yang Li <yang.li@rug.nl>
# Version: 1.0.0
# Date: 10 Dec. 2007


conditionCombination <- function( nEnvFactors, nLevels, Level, envFactorNames )
{
    #This is use by designScore function and experimentDesignTable fuction
    if( nEnvFactors==0 )
    {
        cond.comb <- NULL
    }
    else
    {
       Level2 <- NULL
        for( i in 1:nEnvFactors)
        {
          Level2 [[i]] <- list(NULL)
              if (is.null(Level[[i]]) | !is.numeric(Level[[i]]))
              {
                    Level2[[i]] <- seq(1:nLevels[i])
              } else{
                    Level2[[i]] <- Level[[i]]
              }
        }
        
        if ( nEnvFactors==1 )
            cond.comb <- as.matrix( expand.grid(Level2[[1]]) )
        if( nEnvFactors==2 )
            cond.comb <- expand.grid( Level2[[1]], Level2[[2]] )
        if( nEnvFactors==3 )
            cond.comb <- expand.grid( Level2[[1]], Level2[[2]],Level2[[3]] )

        if (!is.null(envFactorNames))
        {
            colnames(cond.comb)<- as.character(envFactorNames)
        }
    }
    return( cond.comb )
}