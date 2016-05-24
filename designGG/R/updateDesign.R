# Function name: updateDesign
# Author: Yang Li <yang.li@rug.nl>
# Maintainer: Yang Li <yang.li@rug.nl>
# Version: 1.0.0
# Date: 10 Dec. 2007


updateDesign <- function( array.allocation, condition.allocation, nRILs,nSlides,
                          nEnvFactors, nTuple, bTwoColorArray)
{
    if( bTwoColorArray )
    {
        out                          <- arrayUpdate (array.allocation,
                                           condition.allocation, nRILs, nSlides)
        new.array.allocation         <- out[[1]]
        updated.condition.allocation <- out[[2]]

        if( nEnvFactors == 0 )
        {
            new.condition.allocation   <- NULL
        }else{
            new.condition.allocation   <-
                        conditionUpdate (updated.condition.allocation, nTuple,
                                       bTwoColorArray )
        }
    }
    else
    {
      new.array.allocation     <- NULL
      new.condition.allocation <- conditionUpdate (condition.allocation, nTuple,
                                         bTwoColorArray )
    }
    
    return( list(new.array.allocation, new.condition.allocation) )
}
