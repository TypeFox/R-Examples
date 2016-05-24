BinaryCut <- function(iv, dv, nbins=10, 
                      minBin=ceiling(min(table(dv))/50), 
                      woeDelta=0.1,
                      bins=FALSE, debug=FALSE) {
  # converts a numeric input variable into a factor using supervised cutting strategies
  
  # check validity of inputs
  unique.dv <- unique(dv)
  stopifnot(length(iv) == length(dv)) # lengths of inputs must match
  stopifnot(length(unique.dv)==2) # dependent variable must have two classes
  stopifnot(class(iv) %in% c("numeric","integer"))
  stopifnot(nbins >= 2) # must have at least 2 bins
  
  # validate minBin input
  ValidateMinBin <- function(minBin){
    if ((minBin < 0) | (length(minBin)>1) | (minBin>=length(iv))){
      stop("Invalid value for minBin provided, must have length 1 and be between 0 and length(iv): ", length(iv), ". ",
           "Value provided is ", minBin)
    }
  }
  if (class(minBin) %in% c("integer","numeric")){
    # a numeric value was provided
    ValidateMinBin(minBin)
  } else if (class(minBin) == "function") {
    # a function was provided, apply it
    minBinFunc <- minBin
    minBin <- minBinFunc(iv,dv)
    ValidateMinBin(minBin)
  }
  
  # if there are few unique values then reset number of bins
  uniqueValues <- unique(iv)
  uniqueValues <- uniqueValues[!is.na(uniqueValues)]
  nbins.actual <- length(uniqueValues)
  if (nbins.actual < nbins){
    # there are fewer unique values than bins, return the sorted unique values as bin limits
    nbins <- nbins.actual
    binlimits <- sort(uniqueValues)
    # this will merge values in first two bins since we use include.lowest=TRUE
    # add an extra limit below the first unique value
    newLowerLimit <- binlimits[1] - (binlimits[2]-binlimits[1])
    binlimits <- c(newLowerLimit, binlimits)
  } else {
    # there are more unique values than bins
    # execute binning with data from the smaller class
    # try the first unique value of dv, assume it's the smallest
    idx.small <- unique.dv[1] == dv
    min.table.dv <- min(table(dv))
    if (sum(idx.small) == min.table.dv){
      # assumption was correct, first unique value is smallest class
      dvSmall <- unique.dv[1]
    } else if (sum(!idx.small) == min.table.dv) {
      # assumption was incorrect, second unique value is smallest class
      dvSmall <- unique.dv[2]
    } else {
      stop("Count of dv values in each class does not match indices.")
    }
    binlimits <- quantile(iv[dv==dvSmall], probs=seq(0,1,1/nbins), names=FALSE, na.rm=TRUE)
    binlimits <- unique(binlimits) # ensure that bin limits are unique, accounts for ties
  }
  if (debug){cat("Binlimits before validation [", binlimits, "]\n")}
  
  # test for case where there are fewer than 3 binlimits
  if (length(binlimits)<=2){
    # bins are invalid, try binning again without regard to dependent variable
    binlimits <- quantile(iv, probs=seq(0,1,1/nbins), names=FALSE, na.rm=TRUE)
    binlimits <- unique(binlimits) # ensure that bin limits are unique, accounts for ties
  }
  
  # ensure that bin limits are outside min and max values of independent variable
  ivmin <- min(iv, na.rm=TRUE)
  ivmax <- max(iv, na.rm=TRUE)
  # reset min if necessary
  if (ivmin == ivmax) {
    # min and max are identical, there is a single value. Add new min limit below the min.
    binlimits <- c(ivmin-1, ivmax)
  } else if (binlimits[1] > ivmin){
    # The min limit is greater than the min value, replace the min limit
    binlimits[1] <- ivmin
  }
  # reset max if necessary
  if (binlimits[length(binlimits)] < ivmax){
    binlimits[length(binlimits)] <- ivmax
  }
  
  # cut numeric values into factor
  fiv <- cut(iv, breaks=binlimits, include.lowest=TRUE)
  if (debug){cat("Binlimits after validation [", binlimits, "]\n")}
  
  # merge bins if the count of dependent variable values falls below a threshold
  mbbc <- MergeBinsByCount(fiv, iv, dv, nbins, minBin, binlimits, debug)
  binlimits <- mbbc$binlimits
  fiv       <- mbbc$fiv
  if (debug){cat("Binlimits after merge by count [", binlimits, "]\n")}
  
  # merge bins if the difference in the WOE is below a threshold
  woeList <- Woe(fiv, dv)
  #woeDelta <- minWoeFraction * (max(woeList$woe.levels) - min(woeList$woe.levels))
  mbbw <- MergeBinsByWoe(fiv, iv, dv, nbins, woeDelta, binlimits, debug)
  binlimits <- mbbw$binlimits
  fiv       <- mbbw$fiv
  if (debug){cat("Binlimits after merge by WOE [", binlimits, "]\n")}

  
  # catch the case where there is one bin
  if (length(levels(fiv)) == 1 & length(binlimits) == 2) {
    # introduce a new bin boundary between the existing boundaries
    binlimits[3] <- binlimits[2]
    binlimits[2] <- mean(binlimits[c(1,3)])
    # cut again with new binlimits
    fiv <- cut(iv, breaks=binlimits, include.lowest=TRUE)
  }
  
  if (bins){
    # return a list with a factor of levels and a vector of the bin limits / breaks
    return(list(fiv=fiv, breaks=binlimits))
  } else {
    # return a factor
    return(fiv)
  }
}


MergeBinsByWoe <- function(fiv, iv, dv, nbins, woeDelta, binlimits, debug){
  # recursively merge bins until difference in bin Weight of Evidence exceeds threshold
  if (length(binlimits)==2){
    # single bin, return original bins, don't try to merge
    return(list(binlimits=binlimits, fiv=fiv))
  }
  # replace missing values in factor
  fiv <- CleanNaFromFactor(fiv)
  # compute the weight of evidence
  woelist <- Woe(fiv, dv)
  woeDeltaVec <- abs(diff(woelist$woe.levels))
  if (debug) cat("\nMergeBinsByWoe:\n")
  
  # enter loop, keep iterating until difference in woe for adjacent bins exceeds threshold
  while ((min(woeDeltaVec) < woeDelta) & (length(binlimits) > 3)) {
    # find index of first bin with minimum WOE
    imin <- which.min(woeDeltaVec)
    # determine which bin we want to merge
    if (imin == 1){
      # no bin on left, so choose right, remove second boundary
      binlimits <- binlimits[-2]
    } else if (imin == (length(woelist$woe.levels)-1)) {
      # no bin on right, so choose left, remove second-from-last boundary
      binlimits <- binlimits[-(length(binlimits)-1)]
    } else {
      # remove cutpoint corresponding to smallest bins
      binlimits <- binlimits[-(imin+1)]
    }
    if (debug) cat("  WoeDelta=[",woeDeltaVec, "] min=", imin, "newBinLimits=[", binlimits, "]\n")
    # re-cut the bins
    fiv <- cut(iv, breaks=binlimits, include.lowest=TRUE)
    # replace missing values
    fiv <- CleanNaFromFactor(fiv)
    # update WOE
    woelist <- Woe(fiv, dv)
    woeDeltaVec <- abs(diff(woelist$woe.levels))
  }
  
  # return the factor of bins and the bin limits
  return(list(binlimits=binlimits, fiv=fiv))
}


MergeBinsByCount <- function(fiv, iv, dv, nbins, minBin, binlimits, debug){
  # recursively merge bins until none are smaller than the minBin size, 
  # where the counts are broken out by the dependent variable values
  # must have at least 2 bins (3 limits) in order to enter or continue loop
  loopCount <- 1
  if(debug){
    cat("\nMergeBinsByCount: minBin=", minBin, "#####\n")
    print(table(fiv,dv))
  }
  while ( (min(table(fiv,dv)) < minBin) & (length(binlimits) > 3)){
    # remove the boundary between the smallest bin and it's smaller neighbor
    binSizes <- as.vector(table(fiv)) # vector of bin sizes counting all values
    # get index of smallest bin broken out by dependent variable
    # use modulus operator %% to account for min value in 1st column of dv values or second
    iSmallBin <- which.min(table(fiv,dv)) %% length(levels(fiv))
    # fix the case where iSmallBin is zero, this happens when the min is in the last row, 1st column of counts
    if (iSmallBin == 0){
      iSmallBin <- length(levels(fiv))
    }
    if (debug){cat("\n\n\niSmallBin:", iSmallBin, "\n")}
    if (iSmallBin == 1) {
      # no bin on left, so choose right, remove second boundary
      binlimits <- binlimits[-2]
    } else if (iSmallBin == length(binSizes)) {
      # no bin on right, so choose left, remove second-from-last boundary
      binlimits <- binlimits[-(length(binlimits)-1)]
    } else {
      # bins on left and right, choose smaller, ties will pick first / 1
      if (1 == which.min( binSizes[c(iSmallBin-1, iSmallBin+1)] )){
        # the left bin is smaller, remove boundary
        binlimits <- binlimits[-iSmallBin]
      } else {
        # the right bin is smaller, remove boundary
        binlimits <- binlimits[-(iSmallBin+1)]
      }
    }
    # re-cut the bins
    fiv <- cut(iv, breaks=binlimits, include.lowest=TRUE)
    
    if(debug){print(table(fiv,dv))}
    
    # check for infinite loop
    loopCount <- loopCount + 1 # increment counter
    if (loopCount > nbins){
      warning("Maximum iterations in BinaryCut bin merging.")
      break # exit the loop
    }
  }
  # return the factor of bins and the bin limits
  return(list(binlimits=binlimits, fiv=fiv))
}