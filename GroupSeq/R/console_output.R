"outputBounds" <- function(n, alpha, t, lowerBounds, upperBounds,
                           probDifference, probExit, symmetricBoundsYesNo,
                           spendingFunctionUsed) {
    cat("", "\n")
    cat("", "\n")
    cat("#################################", "\n")
    cat("#                               #", "\n")
    cat("# Output of the computed Bounds #", "\n")
    cat("#                               #", "\n")
    cat("#################################", "\n")
    cat("", "\n")
    cat("", "\n")
    cat("number of interim analyses =",n,"\n")

    ##ouput alpha
    if(symmetricBoundsYesNo==1)
    {
      cat("alpha =",alpha[1],"\n")
    }
    else
    {
      cat("Upper alpha = ",alpha[1],"\n")
      cat("Lower alpha = ",alpha[2],"\n")
    }

    ##output names of spending functions that were used
    FunctionNames <- c("O'Brien-Fleming Type","Pocock Type","Power family: alpha* t^phi",
                       "Hwang-Shih-DeCani fammily","Exact Pocock Bounds")
    if(symmetricBoundsYesNo==1)
    {
      cat(FunctionNames[spendingFunctionUsed[1]],"was used as spending Function.","\n")
      cat("","\n")
    }
    else
    {
      cat("Spending Function for UPPER Bound:",FunctionNames[spendingFunctionUsed[1]],"\n")
      cat("Spending Function for LOWER Bound:",FunctionNames[spendingFunctionUsed[2]],"\n")
      cat("","\n")
    }

    ##output the bounds
    times <- data.frame(" *Times*"=t,check.names=FALSE)
    bounds <- data.frame(" *Lower Bounds*"=lowerBounds," *Upper Bounds*"=upperBounds,check.names=FALSE)
    currentAlpha <- data.frame(" *alpha[i]-alpha[i-1]*"=probDifference, check.names=FALSE)
    cumulativeAlpha <- data.frame(" *cumulative alpha*"=probExit, check.names=FALSE)

    print(cbind(format(times,digits=3),format(bounds,digits=5), format(currentAlpha,digits=5), format(cumulativeAlpha,digits=5) ))

}#end <--*function(...)*


"outputDriftWithBounds" <- function(n, probTotal, drift, expectedStoppingTime,
                                    secondTimeScaleIsUsed, t, t2, t2max,
                                    lowerBounds, upperBounds, probStopping,
                                    probExceedingUpper, probExceedingLower,
                                    confidenceLevel) {
  ###INITIALiZE VARIABLES###
  cumulativeExitProb<-0 #cumulative exit probability
  resultExitProb<-0 #exit probability


  ##compute exit probability and cumulative exit probability
  for(i in 1:n)
  {
    resultExitProb[i] <- probExceedingUpper[i]+probExceedingLower[i]

    if(i==1)
    {
      cumulativeExitProb[i] <- resultExitProb[i]
    }
    else
    {
      cumulativeExitProb[i] <- cumulativeExitProb[i-1] + resultExitProb[i]
    }
  }



  ##Start output##

  cat("\n")
  cat("\n")
  cat("############################", "\n")
  cat("#                          #", "\n")
  cat("# Output Drift with Bounds #", "\n")
  cat("#                          #", "\n")
  cat("############################", "\n")
  cat("\n")
  cat("\n")
  cat("*-----*","\n")
  cat("*Drift:",format(drift,digits=5),"\n")
  cat("*-----*","\n")
  cat("Maximum Information:",t2max,"\n")
  cat("Power is ",confidenceLevel,"\n")
  cat("\n")
  cat("\n")

  ##if no second information time is used output without information time
  if(!secondTimeScaleIsUsed || t2max<=0)
  {
    times <- data.frame(" Times"=t,check.names=FALSE)
    bounds <- data.frame(" Lower Bounds"=lowerBounds," Upper Bounds"=upperBounds,check.names=FALSE)
    exitProbability <- data.frame(" exit Prob"=resultExitProb, check.names=FALSE)
    cumulativeExitProbability <- data.frame(" cumulative Exit Prob"=cumulativeExitProb, check.names=FALSE)

    print( cbind( format(times,digits=3), format(bounds,digits=5),
                       format(exitProbability,digits=5), format(cumulativeExitProbability,digits=5) ) )
    cat("\n")

  }

  else
  {
    ##Output with information time.
    times <- data.frame(" Times"=t,check.names=FALSE)
    info <- data.frame(" Info" =t2,check.names=FALSE)
    bounds <- data.frame(" Lower Bounds"=lowerBounds," Upper Bounds"=upperBounds,check.names=FALSE)
    exitProbability <- data.frame(" exit Prob"=resultExitProb, check.names=FALSE)
    cumulativeExitProbability <- data.frame(" cumulative Exit Prob"=cumulativeExitProb, check.names=FALSE)

    print( cbind( format(times,digits=3), format(info,digits=3), format(bounds,digits=5),
                       format(exitProbability,digits=5), format(cumulativeExitProbability,digits=5) ) )
    cat("\n")
  }

}#end <--*function(...)*


"outputForDifferentDrifts" <- function(n, probTotal, drift,
                                       expectedStoppingTime,
                                       secondTimeScaleIsUsed, t, t2, t2max,
                                       lowerBounds, upperBounds, probStopping,
                                       probExceedingUpper, probExceedingLower) {
  ###INITIALiZE VARIABLES###
  cumulativeExitProb<-0 #cumulative exit probability
  resultExitProb<-0 #exit probability


  ##compute exit probability and cumulative exit probability
  for(i in 1:n)
  {
    resultExitProb[i] <- probExceedingUpper[i]+probExceedingLower[i]

    if(i==1)
    {
      cumulativeExitProb[i] <- resultExitProb[i]
    }
    else
    {
      cumulativeExitProb[i] <- cumulativeExitProb[i-1] + resultExitProb[i]
    }
  }



  ##Start output##
  cat("\n")
  cat("n=",n,",  Drift=",drift,"\n")
  cat("Maximum Information:",t2max,"\n")
  cat("-----------------------------------------------------------------------------------------------\n")

  ##if no second information time is used output without information time
  if(!secondTimeScaleIsUsed || t2max<=0)
  {
    times <- data.frame(" Times"=t,check.names=FALSE)
    bounds <- data.frame(" Lower Bounds"=lowerBounds," Upper Bounds"=upperBounds,check.names=FALSE)
    exitProbability <- data.frame(" exit Prob"=resultExitProb, check.names=FALSE)
    cumulativeExitProbability <- data.frame(" cumulative Exit Prob"=cumulativeExitProb, check.names=FALSE)

    print( cbind( format(times,digits=3), format(bounds,digits=5),
                       format(exitProbability,digits=5), format(cumulativeExitProbability,digits=5) ) )
  cat("-----------------------------------------------------------------------------------------------\n")
    cat("\n")

  }

  else
  {
    ##Output with information time.
    times <- data.frame(" Times"=t,check.names=FALSE)
    info <- data.frame(" Info" =t2,check.names=FALSE)
    bounds <- data.frame(" Lower Bounds"=lowerBounds," Upper Bounds"=upperBounds,check.names=FALSE)
    exitProbability <- data.frame(" exit Prob"=resultExitProb, check.names=FALSE)
    cumulativeExitProbability <- data.frame(" cumulative Exit Prob"=cumulativeExitProb, check.names=FALSE)

    print( cbind( format(times,digits=3), format(info,digits=3), format(bounds,digits=5),
                       format(exitProbability,digits=5), format(cumulativeExitProbability,digits=5) ) )
  cat("-----------------------------------------------------------------------------------------------\n")
    cat("\n")
  }

  #Return cumulativeExitProb[n] - needed for the power in the graph,if Task 3 was chosen
  return(cumulativeExitProb[n])

}#end <--*function(...)*

