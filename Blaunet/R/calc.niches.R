calc.niches <-
function(blauObj, dev.range) {
  
  #initialize data objects
  topbounds <- matrix(0, ncol = ncol(blauObj$dimensions), nrow = ncol(blauObj$memberships))
  lowbounds <- matrix(0, ncol = ncol(blauObj$dimensions), nrow = ncol(blauObj$memberships)) 
  means <- matrix(0, ncol = ncol(blauObj$dimensions), nrow = ncol(blauObj$memberships))

  #calculate top and low boundaries
  for (memCyc in 1:ncol(blauObj$memberships)) {
    for (dimCyc in 1:ncol(blauObj$dimensions)) {

      memRows <- which(blauObj$memberships[,memCyc] == 1)
      dimRows <- blauObj$dimensions[memRows, dimCyc]
      memRows <- memRows[!is.na(dimRows)] #gets rid of the missing values in the relevant dimension
      meanData <- blauObj$dimensions[memRows, dimCyc] #rows for relevant dimension

      if (length(meanData) == 0){ #for when there is no information
        #this can happen in one of two cases:
        #1) no members in the group
        #2) all members of group have NA along the relevant dimension
        means[memCyc, dimCyc] <- NA
        topbounds[memCyc, dimCyc] <- NA
        lowbounds[memCyc, dimCyc] <- NA
      }

      else if (length(meanData) == 1){
        #impute our only information if there's 1 obs for the dimension
        means[memCyc,dimCyc] <- meanData #should be just a number
        topbounds[memCyc,dimCyc] <- meanData
        lowbounds[memCyc,dimCyc] <- meanData
      }

      else if (length(meanData) > 1) {
        meanWeights <- blauObj$weights[memRows,]
        means[memCyc,dimCyc] <- sum(meanData*meanWeights)/sum(meanWeights) 
        # Calculate the standard deviation
        # Information on weighted Standard Deviation found at
        # http://www.sosmath.com/CBB/viewtopic.php?t=2656
        sdDenominator <- ((length(meanWeights) - 1) * sum(meanWeights))/(length(meanWeights))
        sdNumerator <- 0
        for (dataCyc in 1:length(meanData)){
          sdNumerator <- sdNumerator + (meanWeights[dataCyc] * (meanData[dataCyc] - means[memCyc,dimCyc])^2 )
        }
        stdDev <- sqrt(sdNumerator/sdDenominator)
        topbounds[memCyc, dimCyc] <- means[memCyc, dimCyc] + stdDev * dev.range
        lowbounds[memCyc, dimCyc] <- means[memCyc, dimCyc] - stdDev * dev.range
        if (lowbounds[memCyc, dimCyc]<0 & min(dimRows,na.rm=T)>=0) lowbounds[memCyc, dimCyc] <- 0
      }
    }
  }
  blauObj$topbounds <- topbounds
  blauObj$lowbounds <- lowbounds
  
  colnames(blauObj$topbounds) <- colnames(blauObj$dimensions)
  rownames(blauObj$topbounds) <- colnames(blauObj$memberships)
  colnames(blauObj$lowbounds) <- colnames(blauObj$dimensions)
  rownames(blauObj$lowbounds) <- colnames(blauObj$memberships)
  

  #calculate if each node is in a given niche
  blauObj$isInNiche <- matrix(0, nrow = nrow(blauObj$memberships), ncol = ncol(blauObj$memberships))

  #the inside 'apply' takes each row in dimensions and checks if it's within the boundaries
  #the outside 'apply' checks if all elements of each row in the matrix are true
  for (memCyc in 1:nrow(blauObj$lowbounds)){
    blauObj$isInNiche[,memCyc] <- apply(t(apply(blauObj$dimensions, 1, function(x) x >= blauObj$lowbounds[memCyc,] & x <= blauObj$topbounds[memCyc,])), 1, all)
  }

  #overwrite NAs with zeroes
  blauObj$isInNiche[is.na(blauObj$isInNiche)] <- 0

  colnames(blauObj$isInNiche) <- vapply(colnames(blauObj$memberships), function(x) paste(x, "niche", sep="_"), "a")
  rownames(blauObj$isInNiche) <- rownames(blauObj$memberships)
  
  return(blauObj)
}
