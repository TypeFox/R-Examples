# Function name: temperatureStep
# Author: Yang Li <yang.li@rug.nl>
# Maintainer: Yang Li <yang.li@rug.nl>
# Version: 1.0.0
# Date: 10 Dec. 2007


temperatureStep <- function( startTemp, maxTempStep, endTemp, nIterations )
{
    #This is used by designGG main function (Simulated annealing process)
    if( startTemp * maxTempStep^nIterations < endTemp )
    {
        temperature.step <- exp((log(endTemp) - log(startTemp))/nIterations)
    }
    else
    {
        temperature.step <- maxTempStep
    }
    return( temperature.step )
}