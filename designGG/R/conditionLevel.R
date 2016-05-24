# Function name: conditionLevel
# Author: Yang Li <yang.li@rug.nl>
# Maintainer: Yang Li <yang.li@rug.nl>
# Version: 1.0.0
# Date: 10 Dec. 2007

conditionLevel <- function( array.allocation, condition.allocation,
                            condition.combination, nEnvFactors )
{
    #This is used in designScore function
    cond.level <- NULL    
    if( is.null(array.allocation) ) 
    {
        for( i.cond in 1:nrow(condition.allocation) )
        {
            ril       <- which(condition.allocation[i.cond,]!=0)
            for( k in 1:length(ril) )
            {
                cond.level <- rbind(cond.level, condition.combination[i.cond,])
            }
        }
        cond.level <-  as.matrix( cond.level)
    }
    else 
    {
        which.cond  <- apply(condition.allocation, 2, function(x) which(x!=0))
        cond.level  <- as.matrix(condition.combination[which.cond,])

    }

    colnames(cond.level) <- paste("F",seq(1:nEnvFactors),sep="")

    
    return( cond.level )

}