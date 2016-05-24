"integrateByTrapez" <-
function(functionValues, numberOfIntegrationIntervalls, gridSize)
{
  ###Initialize variables###
  sumUp<-0
  area<-0


  ##sum up all function values
  ##All but the first and last function values appear twice in the summation
  sumUp <- sumUp + functionValues[1]
  sumUp <- 2*sum(functionValues[2:numberOfIntegrationIntervalls])
  sumUp <- sumUp + functionValues[numberOfIntegrationIntervalls+1]

  ##Multiply the sum by gridSize/2
  area <- (gridSize/2)*sumUp

  return(area)

}#end <--*function(...)*
