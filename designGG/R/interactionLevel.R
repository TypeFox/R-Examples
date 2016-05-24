# Function name: interactionLevel 
# Author: Yang Li <yang.li@rug.nl>
# Maintainer: Yang Li <yang.li@rug.nl>
# Version: 1.0.0
# Date: 10 Dec. 2007


interactionLevel <- function( genotype.level, condition.level, markerIndex,
                              nEnvFactors )
{
    #This is used by designScore function
    genotype.level <- as.matrix(genotype.level) 
    if( nEnvFactors >= 1 )
    {
        x.QF <- colname <- NULL        
        if( nEnvFactors==1 )
            condition.level <- as.matrix(condition.level)            
        for( i.fac in 1:nEnvFactors )
        {
            x.QF    <- cbind( x.QF,
                         genotype.level[markerIndex,]*condition.level[,i.fac])
            colname <- c(colname, paste("QF", i.fac, sep=""))

        }
        colnames(x.QF) <- colname
        interact.mat   <- x.QF
        
        if( nEnvFactors >= 2 )
        {
            x.FF <- colname <- NULL
            for( i.fac in 1:( nEnvFactors-1) )
            {
                for( j.fac in (i.fac+1):nEnvFactors )
                {
                    x.FF    <- cbind( x.FF,
                               condition.level[,i.fac]*condition.level[,j.fac])
                    colname <- c( colname, paste("F",i.fac,"F",j.fac,sep="") )
                }
            }
            colnames(x.FF)  <- colname
            x.QFF <- colname <- NULL
            for( i.fac in 1:(nEnvFactors-1) )
            {
                for( j.fac in (i.fac+1):nEnvFactors )
                {
                    x.QFF   <- cbind(x.QFF, (genotype.level[markerIndex,]
                                             *condition.level[,i.fac]
                                             *condition.level[,j.fac]))
                    colname <- c(colname,paste("QF",i.fac,"F",j.fac,sep=""))
                }
            }
            colnames( x.QFF ) <- colname
            interact.mat    <- cbind( interact.mat, x.FF, x.QFF )

            if ( nEnvFactors == 3 )
            {
                x.FFF                 <-(condition.level[,1]*condition.level[,2]
                                            *condition.level[,3])
                x.QFFF                <- (genotype.level[markerIndex,]
                                     * condition.level[,1] * condition.level[,2]
                                          *condition.level[,3])
                colnames0             <- colnames(interact.mat)
                interact.mat          <- cbind(interact.mat,x.FFF,x.QFFF)
                colnames(interact.mat)<- c(colnames0,"F1F2F3","QF1F2F3")
            }
        }
    }
    return( interact.mat )
}
