# Function name: arrayUpdate
# Author: Yang Li <yang.li@rug.nl>
# Maintainer: Yang Li <yang.li@rug.nl>
# Version: 1.0.0
# Date: 10 Dec. 2007


conditionAllocation <- function ( selectedRILs, genotype, nConditions, nSlides,
                                     nTuple )
{   
    #This is used in initialDesign function
    nSelectedRILs               <- length(selectedRILs)
    cond.allocation             <- matrix(0, nConditions, nSelectedRILs)
    names.cond                  <- paste( "Condition", 1:nConditions, sep="" )
    names.ind                   <- names(selectedRILs)
    dimnames(cond.allocation)   <- list(names.cond, names.ind)

    genotype                    <- as.matrix(genotype)
    genotype.selected           <- genotype[,selectedRILs]

    for( i in 1:nConditions )
    {
        selected  <- sample(nSelectedRILs,1)   
        while ( sum(cond.allocation[,selected]) > 0 )
        {
            selected   <- sample(nSelectedRILs, 1)   
        }
        cond.allocation[i, selected] <- 1
        loc.temp    <- genotype.selected-genotype.selected[,selected]
        distance    <- apply( abs(loc.temp), 2, sum )
        paired      <- apply( abs(cond.allocation), 2, sum )
        largeDis    <- sort(distance[paired==0],
                                        decreasing=TRUE) [1:(floor(nTuple)-1)]
        for( j in 1:(floor(nTuple)-1) )
        {
            candidate <- which(paired==0 & distance==largeDis[j])
            whichone  <- sample(length(candidate),1)
            cond.allocation[i, candidate[whichone]] <- 1
            paired    <- apply(abs(cond.allocation), 2, sum)
        }
    }

    n.extra <- 2*nSlides-floor(nTuple)*nConditions
    if ( n.extra != 0 )
    {
        for( i.extra in 1:(n.extra) )
        {
            paired              <- apply(cond.allocation, 2, sum)
            candidate           <- which(paired==0)
            whichone            <- sample(length(candidate),1)

            ril.number.perCond  <- apply( cond.allocation,1,
                                            function(x) (length(which(x!=0))) )
            candidate.cond      <- which( ril.number.perCond==floor(nTuple))
            whichcond           <- sample(length(candidate.cond),1)

            cond.allocation[candidate.cond[whichcond],candidate[whichone]] <- 1
        }
    }
   return( cond.allocation )
}