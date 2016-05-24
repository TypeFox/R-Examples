"calculateEqualBounds" <-
function(targetAlpha,provisionallyBounds,n,t2)
{
  ## Initialize variables
  newtonFailed <- FALSE
  accuracy<-10^-7
  epsilon<-0
  noError<-FALSE

  #evaluate start value
  meanOfBounds <- sum(provisionallyBounds)/n

    ## set start values for iteration and lower bounds to -infinity
    startingBounds <- seq(meanOfBounds,meanOfBounds,length=n)
    lowerBounds <- seq(-8,-8, length=n)


    ##---------------------------------------------------------------------##
    ##-------- use the "Sekanten-Verfahren" based on Newton Iteration -----##
    ##---------------------------------------------------------------------##
    ##   ##
    ## calculation obeys following pattern whereby x    converges against  ##
    ## the value we are searching for               k+1           ##

    ##             x  _  x                                                 ##
    ##              k     k-1                                              ##
    ##  x    =  ------------------ * f(x )                                 ##
    ##   k+1     f(x ) - f(x   )        k      ##
    ##              k       k-1                ##


    ##set x    and f(x   )
    ##     k-1        k-1
    ##--------------------
    xkMinusOne <- startingBounds[1]

    ## first calculate fxkMinusOne... ##
    vectorOfResults <- computeAlphaLevel(n,0,t2,lowerBounds,startingBounds, 0,     25)
    ##                                                                      |      ||
    ##                                                                   drift  maximum number of interim analysis
    qpos <- vectorOfResults[[2]]
    qneg <- vectorOfResults[[3]]

    ##...then set it
    fxkMinusOne <- abs( sum(qpos+qneg) - targetAlpha )



    ## choose start point for iteration
    ## the more interim analysis we have the closer we choose our startpoint to xkMinusOne
    ## which is done by xkMinusOne-epsilon
    #  i would be fine with better choices if the Reader knows how to choose them

    ## if we got very many interim analysis we have to limit our epsilon to 10^-10
    if(n > 10)
    {
      epsilon <- 10^(-10)
    }
    else
    {
      epsilon <- 10^(-n)
    }



    ##set x  and f(x )
    ##     k        k
    ##----------------
    xk <- xkMinusOne-epsilon

    ## first calculate fxk... ##
    ##for calculating fxk we will need vectors of equal bounds
    xkUpperBoundsVector <- seq(xk,xk,length=n)
    vectorOfResults <- computeAlphaLevel(n,0,t2,lowerBounds,xkUpperBoundsVector,0,25)
    qpos <- vectorOfResults[[2]]
    qneg <- vectorOfResults[[3]]

    ##...then set it
    fxk <- abs( sum(qpos+qneg) - targetAlpha )



    ##initialize number of loops with purpose of controlling convergence
    numberOfLoops<-0

    ## We do 20 iterations maximally -
    ## if we do not have finished then, we won't have convergence at all
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
      ##for calculating fxkPlusOne we will need a vector of equal bounds
      xkPlusOneUpperBoundsVector <- seq(xkPlusOne,xkPlusOne,length=n)
      vectorOfResults <- computeAlphaLevel(n,0,t2,lowerBounds,xkPlusOneUpperBoundsVector,0,25)
      qpos <- vectorOfResults[[2]]
      qneg <- vectorOfResults[[3]]

      ##...then set it
      fxkPlusOne <- abs( sum(qpos+qneg) - targetAlpha )

      ## check if we reached accuracy
      if ( fxkPlusOne <= accuracy)
      {
        ##accuracy is fulfilled - return xkPlusOneBoundVector
        upperBounds <- xkPlusOneUpperBoundsVector
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

    ##If all 20 loops were made, something went wrong
    if (numberOfLoops==20)
    {
      cat("!!! Convergence Problem while computing exact Pocock Bounds !!!","\n")
      cat("    Computation has been stopped!","\n")
      noError<-FALSE
    }

    if(noError)
    {
      ##return Bounds
      return(upperBounds)
    }
    else
    {
      cat("!!! Error while computing exact Pocock Bounds !!!","\n")
      return(noError)
    }


}#end <--*function*
