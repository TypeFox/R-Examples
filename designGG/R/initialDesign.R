# Function name: initialDesign  
# Author: Yang Li <yang.li@rug.nl>
# Maintainer: Yang Li <yang.li@rug.nl>
# Version: 1.0.0
# Date: 10 Dec. 2007


initialDesign <- function( genotype, nRILs, nSlides, nConditions, nTuple,
                            bTwoColorArray )
{ 
    #This is used by designGG function
    name.ind <- colnames( genotype )

    genotype[is.na(genotype)] <- 0.5
    
    if( bTwoColorArray ) 
    {
        array.allocation <- matrix( 0, nrow = nSlides, ncol = nRILs )
        name.array       <- paste( "Array", 1:nSlides, sep="" )
        dimnames(array.allocation) <- list( name.array, name.ind )
        for( i in 1:nSlides )
        {
            if ( i < (0.5*nRILs)+1 )
            {
                red      <- sample(1:nRILs,1)
                while( sum(abs(array.allocation[,red])) > 0 )
                {
                        red <- sample(1:nRILs, 1 )
                }
                array.allocation[i, red] <- 1
                loc.temp  <- genotype - genotype[,red]
                distance  <- apply( abs(loc.temp), 2, sum )
                paired    <- apply( abs(array.allocation), 2, sum )
                if ( i <= (0.5*nRILs) )
                {
                    candidate <- which( paired==0 &
                                        distance==max(distance[paired==0]) )
                }
                else
                {
                    candidate <- which( distance==max(distance) )
                }
                whichone  <- sample( length(candidate),1 )
                array.allocation[i, candidate[whichone]] <- -1
            }
            else
            {
                red         <- sample(1:nRILs, 1)
                array.allocation[i,red] <- 1
                loc.temp    <- genotype-genotype[,red]
                distance    <- apply(abs(loc.temp), 2, sum)
                candidate   <- which(distance==max(distance))
                whichone    <- sample(length(candidate), 1)
                array.allocation[i, candidate[whichone]] <- -1
            }
        }
        if( nConditions == 1) 
        {
            condition.allocation <- NULL

        } else {
            selectedRILs <- NULL
            for (i.slide in 1:nSlides) {
                selectedRILs <- c(selectedRILs,
                                    which(array.allocation[i.slide,]==1)
                                    ,which(array.allocation[i.slide,]==-1))
            }

            condition.allocation <- conditionAllocation ( selectedRILs, genotype
                                                 ,nConditions, nSlides, nTuple)
        }
    }
    else 
    {
        condition.allocation <- matrix( 0, nrow = nConditions, ncol = nRILs )
        name.cond            <- paste( "Condition", 1:nConditions, sep="" )
        dimnames(condition.allocation) <- list(name.cond, name.ind)
        for( i in 1:nConditions )
        {
            if( (i*floor(nTuple)) < nRILs )
            {
                selected       <- sample(1:nRILs, 1)
                while ( sum(abs(condition.allocation[,selected]))>0 )
                {
                    selected   <- sample(1:nRILs, 1)
                }
                condition.allocation[i, selected] <- 1
                loc.temp    <- genotype-genotype[,selected]
                distance    <- apply( abs(loc.temp), 2, sum )
                paired      <- apply( abs(condition.allocation), 2, sum )
                #largeDis    <- unique( sort(distance[paired==0],
                #                        decreasing=T) )[1:(floor(nTuple)-1)]
                largeDis    <-  sort(distance[paired==0],
                                        decreasing=TRUE) [1:(floor(nTuple)-1)]
           
               
                for( j in 1:(floor(nTuple)-1) )
                {
                    if( (i*floor(nTuple))<nRILs )
                    {
                        candidate <- which(paired==0 & distance==largeDis[j])
                    }
                    else
                    {
                        candidate <- which(distance==largeDis[j])
                    }
                      if(any(condition.allocation[i,candidate]==1)){
                    candidate <-candidate[-which(condition.allocation[i,candidate]==1)]
                    }
                
                    whichone  <- sample(length(candidate),1)
                    condition.allocation[i, candidate[whichone]] <- 1
                    paired    <- apply(abs(condition.allocation), 2, sum)
                }
            }
            else
            {   
                selected    <- sample(1:nRILs, 1)
                condition.allocation[i,selected] <- 1
                loc.temp <- genotype - genotype[,selected]
                distance <- apply( abs(loc.temp), 2, sum )
                #largeDis <- unique( sort(distance,
                #                    decreasing=T))[1:(floor(nTuple)-1)]
                largeDis <-  sort(distance,
                                    decreasing=TRUE)[1:(floor(nTuple)-1)]
                                 
                for( j in 1:(floor(nTuple)-1) )
                {
                    candidate <- which(distance==largeDis[j])
                    if(any(condition.allocation[i,candidate]==1)){
                    candidate <-candidate[-which(condition.allocation[i,candidate]==1)]
                    }
                    whichone  <- sample(length(candidate), 1)
                    condition.allocation[i, candidate[whichone]] <- 1
                }
            }
        }
        n.extra <- nSlides-floor(nTuple)*nConditions
        if ( n.extra != 0 )
        {
          for( i.extra in 1:(n.extra) )
          {

              paired    <- apply(abs(condition.allocation), 2, sum)
              if (all(paired>0))
              {
                i.extra.sample      <- sample( nRILs, 1 )
              }else{
                candidate           <- which(paired ==0)
                i.extra.sample      <- candidate[sample(length(candidate),1)]
              }

              ril.number.perCond  <- apply( condition.allocation,1,
                                            function(x) (length(which(x!=0))) )
              candidate.cond     <- which( condition.allocation[,
                                                    i.extra.sample]==0
                                            & ril.number.perCond==floor(nTuple))
              i.extra.cond        <- candidate.cond[sample(
                                        length(candidate.cond),1 )]
              condition.allocation[i.extra.cond,i.extra.sample] <- 1
            }
        }
   #apply(condition.allocation,1,function(x) length(which(x==1)))
   array.allocation <- NULL
    }    
    return( list(array.allocation, condition.allocation) )
}