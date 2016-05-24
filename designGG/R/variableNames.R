# Function name: variableNames 
# Author: Yang Li <yang.li@rug.nl>
# Maintainer: Yang Li <yang.li@rug.nl>
# Version: 1.0.0
# Date: 10 Dec. 2007


variableNames <- function( nEnvFactors, envFactorNames=NULL )
{
    #This is used in designGG function
    if ( is.null(envFactorNames) )
    {
      envFactorNames <- paste( "F", seq(1:nEnvFactors), sep="")
    }
    
     varNames <- "G" #stands for genotype
     
     if( nEnvFactors>=1 )
     {
        varNames <- c( varNames, envFactorNames, paste("Gx", envFactorNames,sep=""))


        if( nEnvFactors >= 2)
        {
            for( i.fac in 1:(nEnvFactors-1) )
            {
                for( j.fac in (i.fac+1):nEnvFactors )
                {
                    varNames <- c(varNames, paste(envFactorNames[ i.fac], "x",
                            envFactorNames[j.fac], sep=""))
                }
            }
            for( i.fac in 1:(nEnvFactors-1))
            {
                for (j.fac in (i.fac+1):nEnvFactors)
                {
                    varNames <- c(varNames, paste("G", "x",
                            envFactorNames[i.fac],"x",envFactorNames[j.fac],
                             sep=""))
                }
            }

            if( nEnvFactors==3 )
            {
                varNames <- c(varNames, paste(envFactorNames[1],"x",
                envFactorNames[2],"x", envFactorNames[3],sep=""),
                paste("Gx",envFactorNames[1],"x",
                envFactorNames[2],"x", envFactorNames[3],sep=""))
            }
        }
    }
    return( varNames)
}