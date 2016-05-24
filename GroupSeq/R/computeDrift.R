"computeDrift" <-
function(n, t, t2, lowerBounds, upperBounds, target, drift, nMax)
{
  probStopping <-0 #is the probability of reaching ith analysis and stopping.
  probExceedingUpper<-0 # the probability of reaching ith and exceeding upper.
  probExceedingLower<-0 #is the probability of reaching ith and exceeding lower.
  probTotal<-0 #is the total probability of rejecting, sum(1,n)(probExceedingUpper+probExceedingLower)

  ##variables from old Fortran Implementation (they are only used if Newton Iteration fails which usually should not happen)
  lo<-0
  hi<-0
  prev<-0
  gotLo<-0
  gotHi<-0
  tol<-1.0e-6
  del<-0.25



  noError<-FALSE
  ##backup start value 'drift', in case Newton Iteration fails
  saveDrift <-drift

  ##---------------------------------------------------------------------##
  ##---- first use the "Sekanten-Verfahren" based on Newton Iteration ---##
  ##---------------------------------------------------------------------##
  ## ##
  ## calculation obeys following pattern whereby x    converges against  ##
  ## the value we are searching for               k+1 ##

  ##             x  _  x                                                 ##
  ##              k     k-1                                              ##
  ##  x    =  ------------------ * f(x )                                 ##
  ##   k+1     f(x ) - f(x   )        k    ##
  ##              k       k-1              ##

    ## Initialize values
    newtonFailed <- FALSE
    numberOfInterimAnalysis <- n


    ##set x    and f(x   )
    ##     k-1        k-1
    ##--------------------
    xkMinusOne <- drift

    ## first calculate fxkMinusOne... ##
    vectorOfResults <- computeAlphaLevel(n, t, t2, lowerBounds, upperBounds, xkMinusOne, nMax)
    vectorOfProbabilities <- vectorOfResults[[2]]

    ##...then set it
    fxkMinusOne <- abs( sum(vectorOfProbabilities)-target )


    ## Choose start point for iteration -
    ## the more interim analysis we have the closer we choose our startpoint to yb[i-1]
    ## which is done by yb[i-1]-epsilon.

    ## This choice leads to an average number of about 4 iterations to converge -
    ## in case of equally spaced interim analysis it is mostly even better leading to about 3 iterations.

    ## There may be better choices than that.

    ## If we got very many interim analysis we have to limit our epsilon to 10^-10
    if(numberOfInterimAnalysis > 10)
    {
      epsilon <- 10^(-10)
    }
    else
    {
      epsilon <- 10^(-numberOfInterimAnalysis)
    }

    ##set x  and f(x )
    ##     k        k
    ##----------------
    xk <- xkMinusOne-epsilon

    ## first calculate fxk... ##
    vectorOfResults <- computeAlphaLevel(n, t, t2, lowerBounds, upperBounds, xk, nMax)
    vectorOfProbabilities <- vectorOfResults[[2]]

    ##...then set it
    fxk <- abs( sum(vectorOfProbabilities)-target )

    ##number of loops with purpose of controlling convergence
    numberOfLoops<-0

    ## We do 20 iterations at maximum
    ## if we do not have finished then, we won't converge at all
    for (j in 1:20)
    {
      numberOfLoops <- j

      ## get new xkPlusOne like shown above
      xkPlusOne <- xk - ( (xk - xkMinusOne)/(fxk - fxkMinusOne) * fxk )

      ##catch xkPlusOne is NOT defined for any reason for example if (fxk-fxkMinusOne == 0)
      if(is.nan(xkPlusOne))
      {
        newtonFailed <- TRUE
        break
      }

      ##catch diverging xkPlusOne
      if(xkPlusOne==Inf || xkPlusOne==-Inf)
      {
        newtonFailed <- TRUE
        break
      }

      ##calculate new fxkPlusOne
      vectorOfResults <- computeAlphaLevel(n, t, t2, lowerBounds, upperBounds, xkPlusOne, nMax)
      vectorOfProbabilities <- vectorOfResults[[2]]
      fxkPlusOne <- abs( sum(vectorOfProbabilities)-target )

      ## check if we reached tolerance
      if ( fxkPlusOne <= tol)
      {
        ##tolerance is fulfilled - return xkPlusOne
        drift<-xkPlusOne
        noError<-TRUE
        break  #leave the loop
      }
      else
      {
        ## not within tolerance yet - set new values and do next iteration
        xkMinusOne<- xk
        fxkMinusOne <- fxk

        xk <- xkPlusOne
        fxk <- fxkPlusOne
      }
    }#end <--*for (j in 1:20)*

    ##If all 20 loops were made, something must be wrong -
    ##therefore we try old method from Fortran Implementation
    if (numberOfLoops==20)
    {
      newtonFailed <- TRUE
    }

###########################################################################
########### Old Implementation - usually it should NOT be used ############
###########################################################################


  ##"Sekanten-Verfahren" failed - so we try old conservative method
  if(newtonFailed)
  {
    ##--------------------------------------------------------##
    ##--Compute upper exit prob for drift, should be target.--##
    ##--------------------------------------------------------##
    count<-0
    probTotalWithinTolerance<-FALSE
    while(!probTotalWithinTolerance)
    {
      count<-count+1
      vectorOfResults <- computeAlphaLevel(n, t, t2, lowerBounds, upperBounds, drift, nMax)
      probExceedingUpper <- vectorOfResults[[2]]
      probTotal<-0

      ##THIS OUGHT TO BE probTotal + probStopping FOR POWER!  It works for
      ##positive drifts since probExceedingLower is tiny but for drift < 0
      ##this ought to fail miserably.
      probTotal <- sum(probExceedingUpper)

      if( abs(probTotal-target) <= tol )
      {
        ##probTotal is within tolerance of target ->leave while loop
        probTotalWithinTolerance<-TRUE  ##leave while-loop
        noError<-TRUE
      }
      #check convergence
      else if( abs(drift-prev) <= tol/10 )
           {
             ##drift changes by less than tol/10, stop calculating
             cat(" Convergence problem in function 'computeDrift' during computing the drift.","\n")
             cat(" !!!Calculation stopped now!!!","\n")
             break
           }

           else if( probTotal>target+tol)
                {
                  ##Make sure hi gives probTotal<=target+tol
                  hi <- drift
                  drift <- drift-del
                  gotHi<-1
                }

                else
                {
                  ##Make sure lo gives probTotal<=target+tol
                  lo<-drift
                  drift<-drift+del
                  gotLo<-1
                }

      ##If bracketing values have been found, bisect.
      if(gotLo==1 && gotHi==1)
      {
        prev <- drift
        drift <- (lo+hi)/2
      }

    }#end <--*while(!probTotalWithinTolerance)*
  }

  ##if one of the routines above worked correctly -> return the calculated value...
  if (noError)
  {
    ##Return value of drift
    return(drift)
  }

  ##...else abort analysis
  else
  {
    print("Error in function 'computeDrift': not converging. Abort analysis!",quote=FALSE)
    return(noError)
  }





}#end <--*function*
