"tailProbability" <-
function(upperBound, previousDensity, numberOfIntegrationIntervalls,lowerIntegrationLimit,upperIntegrationLimit,standardDeviation)
{
  ###Initialize variables###
  tempValue<-0
  grid<-0
  previousGrid<-0

  ##previous grid size
  previousGrid <- (upperIntegrationLimit-lowerIntegrationLimit)/numberOfIntegrationIntervalls


  ##Compute function grid points##
  tempValue <- seq(lowerIntegrationLimit,upperIntegrationLimit,length=numberOfIntegrationIntervalls+1)
  tempValue <- previousDensity * ( 1 - pnorm( (upperBound-tempValue)/standardDeviation ) )


  #Numerical integration
  result <- integrateByTrapez(tempValue, numberOfIntegrationIntervalls, previousGrid)

  return(result)
}#end <--*function(...)*
