# Function name: designScore
# Author: Yang Li <yang.li@rug.nl>
# Maintainer: Yang Li <yang.li@rug.nl>
# Version: 1.0.0
# Date: 10 Dec. 2007


designScore <- function ( genotype, array.allocation, condition.allocation,
                          nEnvFactors, nLevels, Level, nConditions, weight=1,
                          optimality="A", bTwoColorArray, envFactorNames )
{
    # This is used by designGG main function
    genotype.level      <-  NULL
    if( bTwoColorArray )
    {
        genotype <- as.matrix(genotype)
        for( i in 1:nrow(array.allocation) )
        {
            RIL.red         <-  which(array.allocation[i,] == 1)
            RIL.green       <-  which(array.allocation[i,] == -1)
            genotype.level  <-  cbind(genotype.level,
                                    genotype[,c(RIL.red,RIL.green)])
        }
    }else{
        ril.all <- NULL
        for( i.cond in 1:nrow(condition.allocation) )
        {
            ril             <-  which(condition.allocation[i.cond,]!=0)
            ril.all         <-  c(ril.all,ril)
        }
        genotype.level      <- genotype[,ril.all]
    }
    genotype.level[is.na(genotype.level)] <-  0.5

    ril.names1 <- colnames(genotype.level)
    if( nEnvFactors == 0 )
    {
        genotype.matrix <- pairLevel( t(genotype.level),ril.names1 )
        sc              <-  apply(genotype.matrix, 2,
             function(x) (length(x)+sum(x^2))/(length(x)*sum(x^2)-sum(x)^2) )
        sc.temp         <- sum(sc)
        if(sc.temp==0)   sc.temp <- 1e-20
        SC              <-  1/ sc.temp

    }else{
        condition.combination  <- conditionCombination ( nEnvFactors, nLevels,
                                         Level,envFactorNames)
        condition.level        <- conditionLevel ( array.allocation,
                                    condition.allocation, condition.combination,
                                    nEnvFactors)


        SC <- 0
        for( markerIndex in 1:nrow(genotype.level) )
        {
            interact.level     <-  interactionLevel ( genotype.level,
                                    condition.level, markerIndex, nEnvFactors)
            if( bTwoColorArray )
            {
                ril.names1               <- colnames(condition.allocation)
                genotype.matrix          <- pairLevel(genotype.level[markerIndex
                                            ,], ril.names1)
                env.matrix               <- pairLevel( condition.level,
                                                      ril.names1 )
                colnames(env.matrix)     <- colnames(condition.level)
                interact.matrix          <- pairLevel( interact.level,
                                                        ril.names1)
                colnames(interact.matrix)<-  colnames(interact.level)
            }else{
                genotype.matrix <- t(as.matrix( genotype.level[markerIndex,] ))
                env.matrix      <- as.matrix(condition.level)
                interact.matrix <- as.matrix(interact.level)
            }

            design.matrix              <-  cbind( rep(1,nrow(genotype.matrix)),
                                  genotype.matrix, env.matrix, interact.matrix )
                                                       
            colnames(design.matrix)[1:2] <- c("Intercept","Q")
            rownames(design.matrix)    <- paste(rownames(genotype.matrix),
                                         "onArray",seq(1,nrow(genotype.matrix)),
                                         sep="" )
                                         
            if( optimality == "A" )
            {
                X       <-  as.matrix(design.matrix)
                xtx     <-  t(X) %*% X

                singular <- function (x)
                {
                    return (sum(round(eigen(x, TRUE, TRUE)$values, 5) == 0) > 0)
                }

                if( singular(xtx) )
                {
                    sc  <-  0
                }else{
                    sc.temp  <-  sum( weight*diag(solve(xtx)) )
                    if( sc.temp==0 )  { sc.temp <- 1e-20}
                    sc <- 1/sc.temp
                }
            }
            if( optimality == "D" )
            {
                x       <- as.matrix(design.matrix)
                xtx     <- t(x) %*% x
                sc      <- det(xtx)
            }
            SC <- SC + sc
        }
    }

    return(SC)
}
