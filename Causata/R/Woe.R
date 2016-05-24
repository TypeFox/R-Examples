Woe <- function(iv, dv, ...) {
  # generic function to compute the weight of evidence
  UseMethod("Woe", iv)
}

Woe.factor <- function(iv, dv, maxOdds=10000, civ=NULL, ...){
  
  # check for valid inputs
  stopifnot(class(iv) == "factor") # indepenedent variable is a factor
  stopifnot(length(iv) == length(dv)) # length of inputs must match
  stopifnot(length(unique(dv))== 2) # dependent variable must have 2 classes
  stopifnot(any(!is.na(dv))) # missing values are not allowed in dv
  if (!is.null(civ)){
    stopifnot(length(iv) == length(civ)) # continuous variable same length as factor
    stopifnot(class(civ) %in% c("numeric","integer")) # require that civ is numeric
  }
  
  # replace missing values with a "BLANK" level
  iv <- CleanNaFromFactor(iv)
  
  # return a matrix / table of counts of true values of DV
  m <- as.matrix(table(iv, dv))
  # we assume that the second column is a count of TRUE values for dependent variable
  # find where denominator is zero, replace values with 1 to avoid divide by zero
  id0 <- as.logical(m[, 1] == 0) # denominator is false values, 1st column
  if (sum(id0) > 0){
    m[id0, 1] <- 1
  }
  # compute odds
  odds <- as.vector(m[, 2] / m[, 1]) #  #true/#false
  # replace calculated odds where the denominator was zero
  if (sum(id0) > 0){
    odds[id0] <- maxOdds
  }
  # replace calculated odds with limit values where the limits are exceeded
  imax <- odds > maxOdds
  imin <- odds < 1/maxOdds
  if (sum(imax) > 0){
    odds[imax] <- maxOdds
  }
  if (sum(imin) > 0){
    odds[imin] <- 1/maxOdds
  }
  
  # compute log odds, this vector has a length equal the number of levels in iv
  logodds <- log(odds)
  
  # create a vector of log-odds values so that we have an array of values of equal length to iv
  # loop for each level of factor
  logodds.iv <- rep(0, length(iv))
  level.vec <- levels(iv)
  bin.count  <- rep(0, length(level.vec))
  true.count <- rep(0, length(level.vec))
  for (i in 1:length(level.vec)){
    idx <- level.vec[i] == iv
    if (sum(idx)>0){
      logodds.iv[idx] <- logodds[i]
      bin.count[i] <- sum(idx) # count of data points in this bin
      true.count[i] <- m[i, 2] # count of true dv values in this bin
    }
  }
  names(logodds) <- level.vec
  
  # compute additionl metrics related to information value
  prob.i.true  <- m[, 2] / sum(m[, 2]) # probability that iv is in level given that dv is true
  prob.i.false <- m[, 1] / sum(m[, 1]) # probability that iv is in level given that dv is false
  ldr <- prob.i.true / prob.i.false
  ldr[ldr >   maxOdds] <-   maxOdds
  ldr[ldr < 1/maxOdds] <- 1/maxOdds
  log.density.ratio <- log(ldr)
  information.value <- (prob.i.true - prob.i.false) * log.density.ratio
  
  # if continuous data was provided then compute linearity
  if (is.null(civ)){
    linearity <- NA
    civ.woe <- NA
    civ.bin.median <- NA
    civ.bin.mean <- NA
  } else {
    linearity.list <- Linearity(iv, dv, civ, logodds) 
    linearity <- linearity.list$correlation
    civ.woe <- linearity.list$civ.woe
    civ.bin.median <- linearity.list$civ.bin.median
    civ.bin.mean <- linearity.list$civ.bin.mean
  }
  
  # return the odds in a list
  return(list(
    woe.levels = logodds, 
    woe = logodds.iv,
    odds = odds,
    bin.count = bin.count,
    true.count = true.count,
    log.density.ratio = log.density.ratio,
    information.value = information.value,
    linearity = linearity,
    civ.bin.median = civ.bin.median,
    civ.bin.mean = civ.bin.mean,
    civ.woe = civ.woe,
    prob.i.true  = prob.i.true,
    prob.i.false = prob.i.false))
}


Linearity <- function(fiv, dv, civ, woe){
  # compute the weighted correlation in the independent variable and log odds at each bin
  # validate input data
  stopifnot(
    (length(fiv) == length(dv)) & 
    (length(civ) == length(dv)) & 
    (length(woe) == length(levels(fiv))))
  
  # determine if missing values are present in the continuos data
  if (any(is.na(civ))){
    # missing values found, find corresponding level
    missing.level.name <- unique(as.character(fiv[is.na(civ)]))
    stopifnot(length(missing.level.name)==1) # allow only one level for missing values
  } else {
    # set missing level name to some random characters so they won't match an existing level name
    missing.level.name <- "Random characters: )&!@(&% ibfwbefp7238172-"
  }
  
  # loop for each level in the factor
  ib <- 1 # bin counter
  il <- 1 # level counter
  civ.median <- 0
  civ.mean <- 0
  woe.level <- 0
  bin.count <- 0
  for (level in levels(fiv)){
    # if this is the level with missing values then skip it
    if (level != missing.level.name){
      idx <- fiv == level # index of values in this level
      civ.median[ib] <- median(civ[idx]) # median of values at this level
      civ.mean[  ib] <- mean(  civ[idx])
      woe.level[ib] <- woe[il] # woe value at this level
      bin.count[ib] <- sum(idx) # number of values in bin at this level
      ib <- ib + 1 # increment counter
    }
    il <- il + 1 # increment counter
  }
  # compute the weighted correlation in the independent variable and log odds
  # weights are set by the count of records in each bin
  result <- list(
    correlation = corr( cbind(civ.median, woe.level, bin.count)),
    civ.bin.median = civ.median,
    civ.bin.mean = civ.mean,
    civ.bin.count = bin.count,
    civ.woe = woe.level )
  return(result)
}