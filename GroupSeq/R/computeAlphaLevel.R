"computeAlphaLevel" <-
function(n, t, t2, lowerBounds, upperBounds, drift, nMax)
{
  ######################
  #INITIALIZE VARIABLES#
  ######################
  probStopping <-0 #is the probaility of reaching ith analysis and stopping.
  probExceedingUpper<-0 # the probaility of reaching ith and exceeding upper.
  probExceedingLower<-0 #is the probability of reaching ith and exceeding lower.
  expectedStoppingTime<-0 # not implemented yet
  probTotal <-0 #is the total probability of rejecting, sum(1,n)(probExceedingUpper+probExceedingLower)
  gridSize<-0.05 # the grid size for numerical integration by trapezoidal rule
  lastGrid<-0 # the grid of the joint density at the last analysis.
  maxDimForLastGrid<-5000 # the maximum dimension for lastGrid().
  numberOfIntegrationIntervalls<-0 #"integer" the number of intervals for numerical integration
  lowerIntegrationLimit<-0 # the vector of lower integration limits.
  upperIntegrationLimit<-0 # the vector of upper integration limits.
  standardDeviation<-0 #  the standard deviation of the process increment
  standardDeviationProcess<-0 #the standard deviation of the process.
  noError<-TRUE


  ##Check limit for number of interim analyses (looks).
  if(n>nMax)
  {
    cat(" Number of analyses too large - you got ",n,"\n")
    cat(" but allowed are at most ",nMax,"\n")
    cat("change parameter 'nMax' in function 'input.R' for a larger number.","\n")
    break ##stop calculating
  }

  ###--- Calculate standard deviations of increments and process ---###
  ##function stdDeviations(...) returns (standardDeviation, standardDeviationProcess)
  stdDeviationVector<-stdDeviations( n, t2 )
  standardDeviation <- stdDeviationVector[[1]]
  standardDeviationProcess <- stdDeviationVector[[2]]

  ##Set upper (upperIntegrationLimit) and lower (lowerIntegrationLimit) integration limits.

  #upper integration limits
  upperIntegrationLimit <- (upperBounds*standardDeviationProcess)-(drift*t)

  #lower integration limits
  lowerIntegrationLimit <- (lowerBounds*standardDeviationProcess)-(drift*t)

  #number of intervals for numerical integration
  #numberOfIntegrationIntervalls is a vector of length n!
  numberOfIntegrationIntervalls <- abs( trunc( (upperIntegrationLimit-lowerIntegrationLimit) / (gridSize*standardDeviation) ) + 1)

  ###Check limit for number of grid points.
  if( any(numberOfIntegrationIntervalls>maxDimForLastGrid) )
  {
    cat(" Error in function 'computeAlphaLevel'","\n")
    cat(" Bounds too wide - namely (b-a)/gridSize is greater than ",maxDimForLastGrid,"\n")
    break ##stop calculating
  }


  ##----------------------------------------------------##
  #- Begin loop calculating probabilities for each look -#
  #------------------------------------------------------#
  #- For first look some direct calculation is possible -#
  #--- using normal cdf.  Bounds in standard form are ---#
  #-  adjusted for mean of process at t(1). -------------#
  ##----------------------------------------------------##


  ##probability of reaching ith analysis and stopping.
  probStopping[1] <- 1 - ( pnorm( upperBounds[1]-drift*t[1]/standardDeviation[1] )
                  -pnorm( lowerBounds[1]-drift*t[1]/standardDeviation[1] ) )
  ##probability of reaching ith and exceeding upper.
  probExceedingUpper[1] <- 1 - pnorm( upperBounds[1]-drift*t[1]/standardDeviation[1] )

  ##probability of reaching ith and exceeding lower
  probExceedingLower[1] <- pnorm( lowerBounds[1]-drift*t[1]/standardDeviation[1] )


  ##Accumulate boundary crossing probability.
  probTotal <- probTotal + ( probExceedingUpper[1]+probExceedingLower[1] )

  if(n>=2)
  {
    ##After the first look, numerical integration is necessary.
    for (i in 2:n)
    {
      ##-----------------------------------------##
      # Compute density of process at first look. #
      ##-----------------------------------------##

      if(i==2)
      {
        ##call function 'jointDensity' with parameter==1,
        ##for computing joint Density in 1st analysis - look also function jointDensity
        lastGrid <- jointDensity(1, lowerIntegrationLimit[1], upperIntegrationLimit[1], standardDeviation[1],numberOfIntegrationIntervalls,lastGrid)
        #                    |
        #               parameter==1
      }

      ##----------------------------------------------------------------------##
      ##--Compute probStopping, probExceedingUpper and probExceedingLower in function computeCurrentProbability.--##
      ##----------------------------------------------------------------------##
      vectorOfResults <-       computeCurrentProbability(lastGrid,numberOfIntegrationIntervalls,lowerIntegrationLimit,upperIntegrationLimit,i,gridSize,
      standardDeviation[i])
      probStopping[i] <- vectorOfResults[[1]]
      probExceedingUpper[i] <- vectorOfResults[[2]]
      probExceedingLower[i] <- vectorOfResults[[3]]

      probStopping[i] <- 1 - probStopping[i]

      ##If density will be needed for next step, compute it with parameter i
      if(i!=n)
      {
        lastGrid <- jointDensity(i, lowerIntegrationLimit, upperIntegrationLimit, standardDeviation[i], numberOfIntegrationIntervalls, lastGrid)
      }

      ##Accumulate boundary crossing probability.
      probTotal <- probTotal + ( probExceedingUpper[i]+probExceedingLower[i] )


      ## Expected stopping time.  (Not implemented.)
      ## expectedStoppingTime <- expectedStoppingTime + (probExceedingUpper[i]+probExceedingLower[i])*(1.d0-t [i]) ##

    }#end <--*for (i in 2:n)*
  }#end <--*if(n>=2)*

  ## Compute expected stopping time.  (Not implemented.)
  ## expectedStoppingTime = 1 - expectedStoppingTime

  ##Return the values ('probStopping' and 'expectedStoppingTime' are always equal 0 yet, cause of still missing implementations
  toBeReturned<-list(probAndStop=probStopping, probAndExceedingUpper=probExceedingUpper, probAndExceedingLower=probExceedingLower,
  expectedStoppingTime=expectedStoppingTime,totalTypeOneError=probTotal)
  return(toBeReturned)


}#end <--*function(...)*
