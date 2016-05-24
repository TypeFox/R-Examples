# Function name: conditionUpdate 
# Author: Yang Li <yang.li@rug.nl>
# Maintainer: Yang Li <yang.li@rug.nl>
# Version: 1.0.0
# Date: 10 Dec. 2007


conditionUpdate <-function( condition.allocation, nTuple, bTwoColorArray )
{
  #This function is used by UpdateDesign
    exchangeAll <- function( design )
    {
        which.condition             <- sample(nrow(design), 2)
        temp                        <- design[which.condition[1],]
        design[which.condition[1],] <- design[which.condition[2],]
        design[which.condition[2],] <- temp
        design
    }

    exchangeN <- function( design, n )
    {
        which.condition <- sample(nrow(design), 2)

        for (j in 1:n)
        {
            rils.cond1 <- which(design[which.condition[1],]== 1)
            rils.cond2 <- which(design[which.condition[2],]== 1)
            #cat("rils.cond1",rils.cond1,sep="\t",file="test.txt",append=T)
            #cat("\n","rils.cond2",rils.cond2,sep="\t",file="test.txt",append=T)
            flag       <- 1
            nn         <- 0
            while( flag == 1 & nn <= 20 )
            {
                nn <- nn+1
                i <- sample(rils.cond1,1)
                j <- sample(rils.cond2,1)
                if (colnames(design)[i] != colnames(design)[j] &
                        design[which.condition[2],i] ==0
                        & design[which.condition[1],j]==0)
                    flag <- 0
            }
            if(flag==0) {
                design[which.condition[1],i] <- 0
                design[which.condition[1],j] <- 1
                design[which.condition[2],j] <- 0
                design[which.condition[2],i] <- 1
            }
        }
        design
    }
    oneSampleChange <- function (design,which.cond=NULL)
    {
      paired_RIL  <- apply(design, 2, sum)
        un_RIL      <- which(paired_RIL==0)
        if (is.null(which.cond))   which.cond  <- sample(nrow(design), 1)
        cond.RIL    <- which(design[which.cond,]!=0)
        design[which.cond,] <- 0
        if( length(un_RIL)==0 )
        {
            newRIL  <- sample(ncol(design), 1)
            while (any(newRIL==cond.RIL))
            {
                newRIL  <- sample(ncol(design), 1)
            }
        } else {
            newRIL  <- un_RIL[sample(length(un_RIL), 1)]
        }
        new.cond.RIL <- c(sample(cond.RIL,(length(cond.RIL)-1)), newRIL)
        design[which.cond, new.cond.RIL] <- 1

        return( design)
    }
    removeIdenticalCell<-function(design){

       for(i.cond in 1:(nrow(design)-1))
       {
              for( j.cond in (i.cond+1):nrow(design))
              {
                if (all(design[i.cond,]==design[j.cond,]))
                        design <- oneSampleChange(design,j.cond)
                }
        }
        return(design)
    }
    if (bTwoColorArray)
    {
        if (nTuple < 2)
        {
            P  <- c(0,1,1,1)
        }
        else
        {
            P  <- c(0,0.85,0.95,1)
        }
     }else{
        condition.allocation <- removeIdenticalCell (condition.allocation)

        if (nTuple < 2)
        {
            P  <- c(0.5,1,1,1)
        }
        else
        {
            P  <- c(0.1,0.85,0.95,1)
        }
     }
     pr <- runif(1)
     if( pr <= P[1] )
    {
        new.design <- oneSampleChange(condition.allocation)
    }else if( P[1] < pr & pr <= P[2] )
    {
        new.design <- exchangeN(condition.allocation,1)
    }
    else if ( P[2] < pr & pr <= P[3] )
    {
        new.design <- exchangeN( condition.allocation,
                        sample(2:(floor(nTuple)-1),1))
    }
    else
    {
        new.design <- exchangeAll( condition.allocation )
    }

    return( new.design )
}