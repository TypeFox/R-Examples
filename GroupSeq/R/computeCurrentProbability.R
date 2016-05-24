"computeCurrentProbability" <-
function( lastGrid, numberOfIntegrationIntervalls, lowerIntegrationLimit, upperIntegrationLimit, i , gridSize, standardDeviation )
{
  ######################
  #INITIALIZE VARIABLES#
  ######################

  ## Returned variables:
  probStopping <-0 #is the probaility of reaching ith analysis and stopping
                   #here the value of the integral from lowerIntegrationLimit to upperIntegrationLimit.
  probExceedingUpper<-0 # the probaility of reaching ith and exceeding upper.
                   #here the value of the integral from upperIntegrationLimit to infinity.
  probExceedingLower<-0 #is the probability of reaching ith and exceeding lower
                        #here the value of the integral from lowerIntegrationLimit to minus infinity

  ## Local variables:
  valProbStopping<-0 #vector of integrand values for probStopping
  valProbExceedingUpper<-0 #vector of integrand values for probExceedingUpper
  valProbExceedingLower<-0 #vector of integrand values for probExceedingLower
  grid<-0 # grid is the argument for valProbStopping, valProbExceedingUpper and valProbExceedingLower.

  # numberOfIntegrationIntervalls is the number of steps of size gridSize between lowerIntegrationLimit and upperIntegrationLimit.
  gridFromLastStep<-0 #grid size from last step.

  ##Previous grid size.
  gridFromLastStep <- (upperIntegrationLimit[i-1]-lowerIntegrationLimit[i-1]) / numberOfIntegrationIntervalls[i-1]


  ##Function values to be passed to numerical integration
  ##routine are calculated.
  grid <- seq(lowerIntegrationLimit[i-1],upperIntegrationLimit[i-1],length=numberOfIntegrationIntervalls[i-1]+1)
  valProbExceedingUpper <- ( 1 - pnorm( (upperIntegrationLimit[i]-grid) / standardDeviation ) ) * lastGrid
  valProbExceedingLower <- (     pnorm( (lowerIntegrationLimit[i]-grid) / standardDeviation ) ) * lastGrid

  ##--Calls to numerical integration routine.--##

  ## probStopping depends on valProbStopping
  ## probStopping <- integrateByTrapez(valProbStopping, numberOfIntegrationIntervalls[i-1], gridFromLastStep)
  probStopping <- 0

  probExceedingUpper <- integrateByTrapez(valProbExceedingUpper, numberOfIntegrationIntervalls[i-1], gridFromLastStep)
  probExceedingLower <- integrateByTrapez(valProbExceedingLower, numberOfIntegrationIntervalls[i-1], gridFromLastStep)


  ##Return a vector containing the results probStopping, probExceedingUpper, probExceedingLower
  toBeReturned<-list(probAndStop=probStopping, probAndExceedingUpper=probExceedingUpper, probAndExceedingLower=probExceedingLower)
  return(toBeReturned)


}#end <--*function(...)*
