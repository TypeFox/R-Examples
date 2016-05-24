"calculateTask3" <-
function(n,nMax,t,t2,t2max,t3,drift,equallySpacedTimesInput,secondTimeScaleIsUsed, BoundsSymmetry, alpha, phi, usedFunction,TruncateBoundsInput,
         enterBoundsManually, upperBounds, lowerBounds, taskWindow)
{
  ### did user enter bounds himself? otherwise they must be computed before we can go on with calculating the drift ###
  if(!enterBoundsManually)
  {
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
    }

    ##-- asymmetric bounds --##
    else
    {
      upperBounds <- resultsUpperBounds[[2]]
      lowerBounds <- (-1)*resultsLowerBounds[[2]]
    }


    ## if (5) Pocock Type - the real Pocock Bounds' was chosen -
    ## we have to do some extra calculations
    ## The Spending function gives us an approximately Pocock-Design.
    ## To compute the exact Pocock Bounds we will do according to the following pattern:
    ## (1st)we give the bounds with all bounds are equal. As starting value we are using the
    ##      mean of the bounds computed by our Pocock spending function. I figured out that
    ##      in almost every case this is a quite good approximation so far.
    ## (2nd)we compute the probability according to our equal bounds, as we would do, if user
    ##      had chosen Task-3- at the beginning
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
    }#end <--*if(whatSpendingFunctionIsUsed[1]==5 || whatSpendingFunctionIsUsed[2]==5)*

   }#end <--*if(!enterBoundsManually)*



  ##-----------------------------------------------------------##
  ##--Probabilities from bounds, possibly with non zero drift--##
  ##-----------------------------------------------------------##
  vectorOfResults <- computeAlphaLevel(n,t,t2,lowerBounds,upperBounds,drift,nMax)
  probStopping <- vectorOfResults[[1]]
  probExceedingUpper <- vectorOfResults[[2]]
  probExceedingLower <- vectorOfResults[[3]]
  expectedStoppingTime <- vectorOfResults[[4]]
  probTotal <- vectorOfResults[[5]]

  ##output results from function 'computeAlphaLevel'
  guiOutputTask3(n,probTotal,drift,expectedStoppingTime,secondTimeScaleIsUsed,t,t2,t2max,
                 lowerBounds,upperBounds,probStopping,probExceedingUpper,probExceedingLower,
                 BoundsSymmetry,enterBoundsManually,alpha,phi,usedFunction, taskWindow)




}#end <--*function calculateTask3*
