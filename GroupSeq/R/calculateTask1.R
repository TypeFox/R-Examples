"calculateTask1" <-
function(n,t,t2,equallySpacedTimesInput,secondTimeScaleIsUsedInput, BoundsSymmetry, alpha, phi,
         usedFunction,TruncateBoundsInput, taskWindow)
{
  #Initialize Variables
  lowerBounds<-0 # lowerBounds is the vector of lower standardized boundaries
  upperBounds<-0 # upperBounds is the vector of upper standardized boundaries
  probExit<-0 # probExit is a vector of exit probabilities
  probDifference<-0 # probDifference(i) = probExit(i)-probExit(i-1) where probExit is a vector of exit probabilities

  ##Symmetric bounds - call function computeBounds(...)

    if (!BoundsSymmetry==3)
    {
      usedFunction[2]=usedFunction[1]
      results<- computeBounds(n, 0, alpha[1], phi[1], t, t2, BoundsSymmetry, usedFunction[1], TruncateBoundsInput)
    }

    else
    ## Asymmetric bounds - call function computeBounds(...) twice
    ## first time for upper bounds, second time for lower bounds
    {
      resultsUpperBounds <- computeBounds(n, 0, alpha[1], phi[1], t, t2, 1, usedFunction[1], TruncateBoundsInput)
      resultsLowerBounds <- computeBounds(n, 0, alpha[2], phi[2], t, t2, 1, usedFunction[2], TruncateBoundsInput)
    }

    ## get the values depending on one-sided or two-sided test had been made ##
    ##-- symmetric bounds --##
    if (!BoundsSymmetry==3)
    {
      lowerBounds <- results[[1]]
      upperBounds <- results[[2]]
      probExit <- results[[3]]
      probDifference <- results[[4]]
    }

    ##-- asymmetric bounds --##
    else
    {
      upperBounds <- resultsUpperBounds[[2]]
      lowerBounds <- (-1)*resultsLowerBounds[[2]]
      probExit <- resultsUpperBounds[[3]] + resultsLowerBounds[[3]]
      probDifference <- resultsUpperBounds[[4]] + resultsLowerBounds[[4]]
    }

    ## if (5) Pocock Type - the real Pocock Bounds' was chosen -
    ## we have to do some extra calculations
    ## The Spending function gives us an approximately Pocock-Design.
    ## To compute the exact Pocock Bounds we will do according to the following pattern:
    ## (1st)we give the bounds with all bounds are equal. As starting value we are using the
    ##      mean of the bounds computed by our Pocock spending function. I figured out that
    ##      in almost every case this is a quite good approximation so far.
    ## (2nd)we compute the probability according to our equal bounds, as we would do, if user
    ##      had chosen Task -3- at the beginning
    ## (3rd)we use Newton Iteration to adjust the bounds in every Iteration until we get the appropriate alpha

    if(usedFunction[1]==5 || usedFunction[2]==5)
    {
      ##check for symmetric bounds
      #one-sided
      if(BoundsSymmetry==1)
      {
        upperBounds <- calculateEqualBounds(alpha[1],upperBounds,n,t2)
      }

      ##two-sided symmetric
      else if(BoundsSymmetry==2)
           {
             {
               upperBounds <- calculateEqualBounds(alpha[1]/2,upperBounds,n,t2)
               lowerBounds <- -upperBounds
             }
           }

           else
           ## asymmetric bounds -> maybe we have to calculate 2 times
           {
             ##check where (5) Pocock Type - the real Pocock Bounds was chosen
             if(usedFunction[1]==5)
             {
               upperBounds <- calculateEqualBounds(alpha[1],upperBounds,n,t2)
             }
             if(usedFunction[2]==5)
             {
               lowerBounds <- (-1)*calculateEqualBounds(alpha[2],-lowerBounds,n,t2)
             }
           }

    ##-----------------------------------------------------------##
    ##--Probabilities from bounds, possibly with non zero drift--##
    ##-----------------------------------------------------------##

    vectorOfResults <- computeAlphaLevel(n,t2,t2,lowerBounds,upperBounds,0,25)
    probExceedingUpper <- vectorOfResults[[2]]
    probExceedingLower <- vectorOfResults[[3]]

    ## re-compute exit probability and cumulative exit probability
    probDifference<-0
    probExit<-0

    for(i in 1:n)
    {
      probDifference[i] <- probExceedingUpper[i]+probExceedingLower[i]

      if(i==1)
      {
        probExit[i] <- probDifference[i]
      }
      else
      {
        probExit[i] <- probExit[i-1] + probDifference[i]
      }
    }

  }#end <--*if(whatSpendingFunctionIsUsed[1]==5 || whatSpendingFunctionIsUsed[2]==5)*

    guiOutputTask1(n, alpha, phi, t, lowerBounds, upperBounds, probDifference, probExit,
                   BoundsSymmetry, usedFunction, taskWindow)


}#end <--*function(...)*
