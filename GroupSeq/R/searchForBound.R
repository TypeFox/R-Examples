"searchForBound" <-
function (lastGrid, numberOfIntegrationIntervalls, i, gridSize, probDifference, standardDeviation, lowerIntegrationLimit, upperIntegrationLimit, numberOfInterimAnalysis)
{
  ###Initialize variables###
  noError<-TRUE

  ##Initialize tolerance for result - that implies our accuracy
  tolerance<-1.0e-07


    ##---------------------------------------------------------------------##
    ##---- first use the "Sekanten-Verfahren" based on Newton Iteration ---##
    ##---------------------------------------------------------------------##
    ##   ##
    ## calculation obeys following pattern whereby x    converges against  ##
    ## the value we are searching for               k+1   ##

    ##             x  _  x                                                 ##
    ##              k     k-1                                              ##
    ##  x    =  ------------------ * f(x )                                 ##
    ##   k+1     f(x ) - f(x   )        k      ##
    ##              k       k-1                ##


    ## Initialize values
    newtonFailed <- FALSE

    ##set x    and f(x   )
    ##     k-1        k-1
    ##--------------------
    xkMinusOne <- upperIntegrationLimit[i-1]
    fxkMinusOne <- abs(tailProbability(xkMinusOne, lastGrid, numberOfIntegrationIntervalls[i-1], lowerIntegrationLimit[i-1], upperIntegrationLimit[i-1], standardDeviation)-probDifference)

    ## choose start point for iteration
    ## the more interim analysis we have the closer we choose our startpoint to upperIntegrationLimit[i-1]
    ## which is done by (xkMinusOne-epsilon)

    ## this choice leads to an average number of about 4 iterations to converge -
    ## in case of equally spaced interim analysis it is mostly even better leading to about 3 iterations

    ## in spite of that i would be fine with better choices
    ## if the Reader knows how to choose them

    ## if we got very many interim analysis we have to limit our epsilon
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
    fxk <- abs(tailProbability(xk, lastGrid, numberOfIntegrationIntervalls[i-1], lowerIntegrationLimit[i-1], upperIntegrationLimit[i-1], standardDeviation)-probDifference)

    numberOfLoops<-0

###########################################################################
########################### BEGIN SEARCHING ###############################
###########################################################################

    ## We do 20 iterations maximally -
    ## if we do not have finished then, we won't have convergence at all
    for (j in 1:20)
    {
      numberOfLoops <- j

      ## get new xkPlusOne like shown above
      xkPlusOne <- xk - ( (xk - xkMinusOne)/(fxk - fxkMinusOne) * fxk )

      ##catch xkPlusOne is not defined for any reason for example if (fxk-fxkMinusOne == 0)
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
      fxkPlusOne <- abs(tailProbability(xkPlusOne, lastGrid, numberOfIntegrationIntervalls[i-1], lowerIntegrationLimit[i-1], upperIntegrationLimit[i-1], standardDeviation)-probDifference)


      ## check if we reached tolerance
      if ( fxkPlusOne <= tolerance)
      {
        ##tolerance is fulfilled - return xkPlusOne
        upperIntegrationLimit[i]<-xkPlusOne
        noError<-TRUE
        break ##leave the loop
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
    searchingWithOldMethod <- TRUE

    ##Initialize variables and values
    numberOfLoops<-0

    ##Initialize estimates at previous integration limit.
    uppr <- upperIntegrationLimit[i-1]

    ##Initialize step size.
    del<-10.0

    q <- tailProbability(uppr, lastGrid, numberOfIntegrationIntervalls[i-1], lowerIntegrationLimit[i-1], upperIntegrationLimit[i-1], standardDeviation)

    while(searchingWithOldMethod)
    {
       ##--------------------------------------------------------##
       ##--------------------- q is alright ---------------------##
       ##--------------------------------------------------------##
       ##If q and probDifference are nearly equal, set upperIntegrationLimit[i] and return
       if ( (abs(q-probDifference)) <= tolerance)
       {
         upperIntegrationLimit[i]<-uppr
         noError<-TRUE
         searchingWithOldMethod<-FALSE
       }

       ##--------------------------------------------------------##
       ##------------------- q is too large ---------------------##
       ##--------------------------------------------------------##
       ##Else if q is too large, start increasing uppr by steps.
       else if (q>(probDifference+tolerance))
       {
         ##count for-loops for the purpose of controlling convergence
         numberOfLoops<-0

         ##Reduce step size by factor of 10.
         del <- del/10

         ##Increase uppr by del...
         for (j in 1:50)
         {
           numberOfLoops <- numberOfLoops + 1
           uppr <- uppr + del

           #...and check whether q is near probDifference.
           q <- tailProbability(uppr, lastGrid, numberOfIntegrationIntervalls[i-1], lowerIntegrationLimit[i-1], upperIntegrationLimit[i-1], standardDeviation)

           if (q<=(probDifference+tolerance))
           {
             break ##leave the for-loop and do another *while(stillSearching)*-loop
           }

           #If many iterations do not converge, print warning after each 10 loops.
           if ((j %% 10) == 0)
           {
             print(" Large change in bounds, possible error.",quote=FALSE)
             print(" Time/Point of interim analysis:",quote=FALSE)
             print(i)
           }
         }#end <--*for*
         ##If all 50 loops were made, something must be wrong - abort!
         if (numberOfLoops==50)
         {
           searchingWithOldMethod<-FALSE
           noError<-FALSE
         }

       }#end <--*else if (q>(probDifference+tolerance))*

       ##--------------------------------------------------------##
       ##------------------- q is too small ---------------------##
       ##--------------------------------------------------------##
       ##Else if q is too small, start decreasing uppr by steps.
       else if (q<(probDifference-tolerance))
            {
              ##count for-loops for the purpose of controlling convergence
              numberOfLoops<-0

              ##Reduce step size by factor of 10.
              del <- del/10

              ##Increase uppr by del...
              for (j in 1:80)
              {
                numberOfLoops <- numberOfLoops + 1
                uppr <- uppr - del

                #...and check whether q is near probDifference.
                q <- tailProbability(uppr, lastGrid, numberOfIntegrationIntervalls[i-1], lowerIntegrationLimit[i-1], upperIntegrationLimit[i-1], standardDeviation)

                if (q>=(probDifference-tolerance))
                {
                  break ##leave the for-loop and do another *while(searchingWithOldMethod)*-loop
                }

                #If many iterations do not converge, print warning after each 10 loops.
                if ((j %% 10) == 0)
                {
                  print(" Large change in bounds, possible error.",quote=FALSE)
                  print(" Time/Point of interim analysis:",quote=FALSE)
                  print(i)
                }
              }#end <--*for*
              ##If all 80 loops were made, something must be wrong - abort!
              if (numberOfLoops==80)
              {
                searchingWithOldMethod<-FALSE
                noError<-FALSE
              }

            }#end <--*else if (q>(probDifference+tolerance))*

    }#end <--*while(searchingWithOldMethod)*
  }#end <--*if(newtonFailed)*

     else
     {
       ## "Sekanten-Verfahren" was successful
     }


  ##if one of the routines above worked correctly -> return the calculated value...
  if (noError)
  {
    return(upperIntegrationLimit[i])
  }

  ##...else abort analysis
  else
  {
    print("Error in search: not converging. Abort analysis!",quote=FALSE)
    return(noError)
  }



}#end <--*function(...)
