# Function name: arrayUpdate
# Author: Yang Li <yang.li@rug.nl>
# Maintainer: Yang Li <yang.li@rug.nl>
# Version: 1.0.0
# Date: 10 Dec. 2007


arrayUpdate <- function( array.allocation, condition.allocation, nRILs, nSlides)
{
    #This is used by updateDesign function
    samplesExchange <- function( design, cond.design )
    {
        which.array  <- sample(nrow(design), 2)
        array1       <- c(which(design[which.array[1],]==1),
                                        which(design[which.array[1],]==-1) )
        array2       <- c(which(design[which.array[2],]==1),
                                        which(design[which.array[2],]==-1) )
        uniq.RIL     <- unique(c(array1, array2))
        RIL          <- rep(sample(uniq.RIL, length(uniq.RIL)), length=4)
        design[which.array[1],]        <-   0
        design[which.array[2],]        <-   0
        design[which.array[1], RIL[1]] <-   1
        design[which.array[1], RIL[2]] <-  -1
        design[which.array[2], RIL[3]] <-   1
        design[which.array[2], RIL[4]] <-  -1

        
        if (!is.null(cond.design)){
            a1     <- which.array[1]
            a2     <- which.array[2]
            temp   <- cond.design[, c((a1*2-1):(a1*2), (a2*2-1):(a2*2))]
            
            if(length(uniq.RIL)==4){
                temp2  <- temp[,colnames(design)[RIL]]
            } else {      
                RIL.old               <- c(array1,array2)
                uniq.temp             <- temp[,which(duplicated(RIL.old)==FALSE)]
                duplic.temp           <- as.matrix(temp[,
                                            which(duplicated(RIL.old)==TRUE)])
                uniq.temp2            <- uniq.temp[,colnames(design)
                                          [RIL[which(duplicated(RIL)==FALSE)]]]
                colnames(duplic.temp) <- colnames(design)[
                                          RIL[which(duplicated(RIL)==TRUE)]]
                temp2                 <- cbind(uniq.temp2, duplic.temp)
            }
            cond.design[, c((a1*2-1):(a1*2), (a2*2-1):(a2*2))]         <- temp2
            colnames(cond.design)[c((a1*2-1):(a1*2), (a2*2-1):(a2*2))] <-
                                                                colnames(temp2)
        }         
        return( list(design, cond.design ))
    }

    twoSamplesChange <- function (design, cond.design)
    {
        paired_RIL   <- apply(abs(design), 2, sum)    
        unpaired_RIL <- which(paired_RIL==0)          
        which.array  <- sample(nrow(design), 1)
        array.RIL    <- which(design[which.array,]!=0)
        design[which.array,] <- 0
        if( length(unpaired_RIL <=2) )
        {
            newRIL <- sample(ncol(design), 2)
            while ( any(newRIL == array.RIL) )    
            {
                newRIL <- sample(ncol(design), 2)
            }
        }else{
            newRIL <- sample(unpaired_RIL , 2)
        }
        design[which.array, newRIL[1]] <- 1
        design[which.array, newRIL[2]] <- -1
        
        if (!is.null(cond.design)){
            colnames(cond.design) [(which.array*2-1):(which.array*2)] <-
                                                    colnames(design)[newRIL]
        }
        return( list(design, cond.design ))
    }

    oneSampleChange <- function (design,cond.design)
    {
        paired_RIL  <- apply(abs(design), 2, sum)    
        un_RIL      <- which(paired_RIL==0)
        which.array <- sample(nrow(design), 1)
        array.RIL   <- which(design[which.array,]!=0)
        design[which.array,] <- 0
        if( length(un_RIL)==0 )
        {
            newRIL  <- sample(ncol(design), 1)
            while (any(newRIL==array.RIL))       
            {
                newRIL  <- sample(ncol(design), 1)
            }
        } else {
            newRIL  <- sample(un_RIL, 1)
        }
        new.array.RIL <- sample(c(sample(array.RIL, 1), newRIL), 2)
        design[which.array, new.array.RIL[1]] <- 1
        design[which.array, new.array.RIL[2]] <- -1
        
        if (!is.null(cond.design)){
          colnames(cond.design) [(which.array*2-1):(which.array*2)] <-
                colnames(design)[new.array.RIL]
        }
        return( list(design, cond.design ))
    }

    pairUpdate <- function(design, cond.design, P)
    {
        pr <- runif(1)
        if( pr <= P[1] )
        {
            design  <- twoSamplesChange(design, cond.design)
        }
        else if( P[1]<pr & pr<=P[2] )
        {
            design  <- oneSampleChange(design,cond.design)
        }
        else {
            design  <- samplesExchange(design,cond.design)
        }
        design
    }

    if ( nSlides == ceiling( nRILs*0.5 ) )
    {
        p       <- c(0,0,1)
        subset  <- 0
    }else if( nSlides < ceiling( nRILs*0.5 ) )      
    {
        p       <- c(0.2, 0.8, 1)
        subset  <- 0
    }else{
        p       <- c(0.4, 0.9, 1)
        subset  <- ceiling(nRILs*0.5)
        nSubset <- floor( nSlides /(0.5*nRILs) )
    }

    if ( subset==0 )
    {
        out                      <- pairUpdate(array.allocation,
                                            condition.allocation,p)
        new.array.allocation     <- out[[1]]
        new.condition.allocation <- out[[2]]
    } else if( subset > 0 )
        {
        for( i in 1:nSubset){
            varName <- paste("design", i ,sep="")
            assign( varName, array.allocation[((i-1)*subset+1):(i*subset),] )

            varName2 <- paste("condition.allocation", i, sep="")
            if(!is.null(condition.allocation)){
                assign( varName2,
                        condition.allocation[,(2*(i-1)*subset+1):(2*i*subset)] )
            }  else {
                assign( varName2, NULL )
            }
        }
        if( (nSubset*subset) < nSlides){
          design.left               <- array.allocation[(nSubset*subset+1):nSlides,]
          condition.allocation.left <- NULL                                         #mar 9 2011 YL
          if(!is.null(condition.allocation)){
              condition.allocation.left <- condition.allocation[,(2*nSubset
                                              *subset+1):(2*nSlides)]
              }
        }
        new.array.allocation        <- array.allocation
        new.condition.allocation    <- condition.allocation

        if(sample( nSlides, 1) <= (nSubset * subset ) )
        {
            k                           <- sample(nSubset,1)
            design.temp                 <- get(paste("design",k,sep=""))
            condition.allocation.temp   <- get(paste("condition.allocation", k,
                                                    sep=""))
            out                         <- samplesExchange(design.temp,
                                                    condition.allocation.temp)
            design.temp2                <- out[[1]]
            condition.allocation.temp2  <- out[[2]]

            new.array.allocation[((k-1)*subset+1):(k*subset),] <- design.temp2
            if (!is.null( new.condition.allocation)){
                new.condition.allocation[,(2*(k-1)*subset+1):(2*k*subset)]  <-
                                                    condition.allocation.temp2
            colnames(new.condition.allocation)[(2*(k-1)*subset+1):
                          (2*k*subset)] <- colnames( condition.allocation.temp2)
            }

        }else{
            if ( nrow(design.left) < 2 )  p <- c(0.6, 1, 1)
            out                         <- pairUpdate(design.left,
                                              condition.allocation.left,p)
            design.left2                <- out[[1]]
            new.array.allocation[(nSubset*subset+1):nSlides,] <- design.left2

            condition.allocation.left2  <- out[[2]]
            if (!is.null( new.condition.allocation)){
                new.condition.allocation[,(2*nSubset*subset+1):(2*nSlides)] <-
                                                condition.allocation.left2
                colnames(new.condition.allocation)[(2*nSubset*subset+1):
                        (2*nSlides)] <- colnames(condition.allocation.left2)
            }
        }
    }
    return( list(new.array.allocation, new.condition.allocation) )
}