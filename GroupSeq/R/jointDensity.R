"jointDensity" <-
function (i, lowerIntegrationLimit, upperIntegrationLimit,standardDeviation, numberOfIntegrationIntervalls, lastGrid)
{
  ##first analysis
  if (i==1)
  {
    ###Initialize variables###
    currentGridSize<-0


    ##Compute grid size to be used
    currentGridSize <- (upperIntegrationLimit[1]-lowerIntegrationLimit[1])/numberOfIntegrationIntervalls[1]


    ##Evaluate function (normal density) at grid points
    ##lastGrid <- seq(lowerIntegrationLimit,upperIntegrationLimit,by=currentGridSize)
    lastGrid <- seq(lowerIntegrationLimit,upperIntegrationLimit,length=numberOfIntegrationIntervalls[1]+1)
    lastGrid <- dnorm( lastGrid/standardDeviation ) / standardDeviation

  }

  ##second analysis and later
  else
  {
    ###Initialize variables###
    currentGridSize<-0 # 'currentGridSize' is the grid size used.
    grid<-0 # 'grid' is the argument for function evaluation.
    lastCopy<-0 # 'lastCopy' is a temporary vector of function values because 'lastGrid' from
        # previous step itself is needed to compute each current value.

    ##Current grid size to be used
    currentGridSize <- (upperIntegrationLimit[i]-lowerIntegrationLimit[i])/numberOfIntegrationIntervalls[i]

    ##Previous grid size.
    previousGridSize <- (upperIntegrationLimit[i-1]-lowerIntegrationLimit[i-1]) / numberOfIntegrationIntervalls[i-1]

    lastlast <- seq(lowerIntegrationLimit[i-1],upperIntegrationLimit[i-1],length=numberOfIntegrationIntervalls[i-1]+1)


    ##Evaluate function over grid lowerIntegrationLimit + [j-1]*currentGridSize, j=1,numberOfIntegrationIntervalls+1.
    for ( j in 1 : (numberOfIntegrationIntervalls[i]+1) )
    {
      grid <- lowerIntegrationLimit[i] + ( currentGridSize * (j-1) )
      f <- lastGrid * (dnorm((grid-lastlast)/standardDeviation) / standardDeviation )
      lastCopy[j] <- integrateByTrapez ( f, numberOfIntegrationIntervalls[i-1], previousGridSize )
    }


    ##Return vector 'lastGrid' respectively 'lastCopy'
    return(lastCopy)

  }#end <--*else*

}#end <--*function(...)*
