FragMatch <- function(newData, refData, maxValue = 1000, errorBound = 25,
                      weight = 1, na.rm = TRUE){
  if(length(maxValue) != length(errorBound))
    stop("'maxValue' and 'errorBound' must have identical lengths.")
  if(any(maxValue <= 0))
    stop("'maxValue' has to be a vector of positive values.")
  if(any(errorBound <= 0))
    stop("'errorBound' has to be a vector of positive values.")
  if(length(weight) > 1){
    weight <- weight[1]
    warning("Only first element of 'weight' is used.")
  }
  if(weight < 0)
    stop("'weight' has to be a non-negative real.")
  
  errorBound <- errorBound[order(maxValue)]
  maxValue <- sort(maxValue)
  
  refNames <- unique(refData$Sample)
  newNames <- unique(newData$Sample)
  
  res <- matrix("", nrow = length(newNames), ncol = length(refNames))
  rownames(res) <- newNames
  colnames(res) <- refNames
  
  for(i in 1:length(newNames)){
    for(j in 1:length(refNames)){
      New <- newData[newData$Sample == newNames[i],]      
      Ref <- refData[refData$Sample == refNames[j],]
      nrEnz <- length(unique(Ref$Enzyme))
      res.enz.match <- max.pos.match <- numeric(nrEnz)
      refenzNames <- unique(Ref$Enzyme)
      newenzNames <- unique(New$Enzyme)
      
      if(length(refenzNames) != length(newenzNames))
        stop("Number of enzymes is different for new data and reference data!")
      if(!all(refenzNames == newenzNames))
        stop("Names of enzymes are different for new data and reference data!")
      
      for(k in 1:nrEnz){
        Newk <- New[New$Enzyme == newenzNames[k],]
        Refk <- Ref[Ref$Enzyme == newenzNames[k],]
        ## number of matches
        max.pos.match[k] <- nrow(Refk)
        res.enz.match[k] <- .countMatches(abs(Newk$MW), abs(Refk$MW), 
                                          maxValue = maxValue, 
                                          errorBound = errorBound)
      }
      
      ind.weights <- max.pos.match == res.enz.match
      max.match <- sum(weight*max.pos.match)
      res.match <- sum(ifelse(ind.weights, weight*res.enz.match, res.enz.match))
      res[i,j] <- paste(res.match, "_", max.match, sep = "")
    }
  }
  res
}

.countMatches <- function(x, y, maxValue, errorBound){
  nrMatch <- 0
  if(any(y > max(maxValue)) | any(x > max(maxValue))){
    errorBound <- c(errorBound, errorBound[which.max(maxValue)])
    maxValue <- c(maxValue, max(y,x))
  }
  for(i in 1:length(maxValue)){
    yi <- y[y < maxValue[i]]
    xi <- x[x < (maxValue[i]+errorBound[i])]
    for(j in 1:length(yi)){
      ind <- (xi > (yi[j]-errorBound[i])) & (xi < (yi[j]+errorBound[i]))
      if(any(ind)){
        nrMatch <- nrMatch + 1
        xi <- xi[-which(ind)[1]]
      }
    }
  }
  nrMatch
}
