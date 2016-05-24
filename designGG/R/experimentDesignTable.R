# Function name: experimentDesignTable
# Author: Yang Li <yang.li@rug.nl>
# Maintainer: Yang Li <yang.li@rug.nl>
# Version: 1.0.0
# Date: 10 Dec. 2007

experimentDesignTable <- function( array.allocation, condition.allocation,
                                   nEnvFactors, nLevels, Level, fileName,
                                   envFactorNames, directory )
{
    if( !is.null(array.allocation) )
    {
        if(is.null(colnames(array.allocation)))
        {
            colnames(array.allocation) <- paste("Strain",
                                        seq(1,ncol(array.allocation)),sep="")
        }
        
        arrayDesign <- cbind( colnames(array.allocation)
                            [apply(array.allocation,1,function(x) which(x==1))],
                              colnames(array.allocation)
                           [apply(array.allocation,1,function(x) which(x==-1))])
        colnames(arrayDesign) <- c(colnames(arrayDesign),"Channel1","Channel2")
        arrays     <- paste("array",seq(1,nrow(array.allocation)),
                                        sep="")
        arrayDesign <- cbind(arrays,arrayDesign)
        colnames(arrayDesign)[1] <-  "Array"
        if(!is.null(fileName))
            {write.csv(arrayDesign, file=paste(directory,"/",fileName,"_arrayDesign.csv",
                    sep=""),row.names=FALSE)
            }

    } else  {
        arrayDesign <- NULL
    }

    if(nEnvFactors>0){
        condition.combination <- conditionCombination( nEnvFactors, nLevels,
                                                     Level, envFactorNames)
        for(i in 1:nEnvFactors){
            if( is.character(Level[[i]]) )
                  for(j in 1:nLevels[i]){
                    condition.combination[(condition.combination[,i] ==j),i] <-
                        Level[[i]][j]
                    }
            }
        conditionDesign0          <-  matrix(0,nrow=nrow(condition.combination),
                                       ncol=max(apply(condition.allocation, 1,
                                            function(x) length(which(x!=0)))) )
        colnames(conditionDesign0)<-  paste("StrainName",
                                           seq(1,ncol(conditionDesign0)),sep="")

        if(is.null(colnames(condition.allocation)))
        {
            colnames(condition.allocation) <- paste("Strain",
                                       seq(1,ncol(condition.allocation)),sep="")
        }

        for( i.cond in 1:nrow(condition.allocation) )
        {
            ril     <-   which(condition.allocation[i.cond,]!=0)
            conditionDesign0[i.cond,1:length(ril)]  <-
                                          colnames(condition.allocation)[ ril ]
        }

        conditionDesign0[ which(conditionDesign0==0) ] <- NA
        conditionDesign <- cbind( condition.combination, conditionDesign0 )
        conditions <- paste("condition",
                                     seq(1,nrow(condition.combination)), sep="")
        conditionDesign <- cbind(conditions,conditionDesign)
        colnames(conditionDesign)[1] <- "Condition"
        if(!is.null(fileName)){
            write.csv( conditionDesign, file=paste(directory,"/",fileName,
                                "_conditionDesign.csv", sep=""), row.names=FALSE )
            }
      }  else{
         conditionDesign <- NULL
      }
      
    return ( list(arrayDesign = arrayDesign, conditionDesign = conditionDesign))
}
