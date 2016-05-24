
PredictivePower <- function(iv, ...) {
  # generic function
  UseMethod("PredictivePower", iv)
}


PredictivePower.factor <- function(iv, dv, warn.levels=30, cv=NULL, debug=FALSE, ...){
  # test inputs
  stopifnot(length(iv) == length(dv)) # lengths of inputs must match
  stopifnot(length(unique(dv)) == 2) # dependent variable must have only two classes
  stopifnot(any(!is.na(dv))) # dependent variable cannot have missing values
  if (length(levels(iv)) >= warn.levels){
    warning("Number of levels in factor ", length(levels(iv)), " exceeds threshold ", warn.levels)
  }
  
  # determine if cross validation is used
  if (is.null(cv)) {
    # cross validation is not used, set indices to use all data
    idx1 <- rep(TRUE, length(iv))
    idx2 <- idx1
  } else {
    # the cv argument is not null, verify that it's an array of booleans
    stopifnot(length(cv) == length(iv))
    stopifnot(class(cv) == "logical")
    # set indices to use opposite sets of data
    idx1 <-  cv
    idx2 <- !cv
  }
  
  # build a table of counts for each level in iv (rows) broken out by dv values (columns)
  tab1 <- table(iv[idx1],dv[idx1])
  tab2 <- table(iv[idx2],dv[idx2])
  # compute rate of positive values (2nd level of dv) within each bin
  rowSumsTab1 <- rowSums(tab1)
  rowSumsTab2 <- rowSums(tab2)
  class1avg1 <- tab1[, 2] / rowSumsTab1
  # for each bin of iv compute the percent of all values in each bin
  pctRec2 <- rowSumsTab2 / sum(rowSumsTab2)
  # for each bin of iv compute the percent of all positive dv values in each bin
  pctDv2 <- tab2[, 2] / sum(tab2[, 2])
  # sort by rate of positives using the first CV set
  sortResult1 <- sort(class1avg1, decreasing=TRUE, index.return=TRUE)
  idxBins1 <- sortResult1$ix
  
  # compute the area under the predictive power curve
  Ar2 <- 0 # running sum of area under gains chart
  dvRunningSum2 <- 0 # running total of the percent of positive dv values
  for (i in 1:length(idxBins1)){
    j <- idxBins1[i] # get the index of bins as sorted by average dv value
    # compute the area of the "box" below the curve using the second set of data
    Abox2 <- dvRunningSum2 * pctRec2[j]
    # area under the gains chart for this bin, added to area from prior bins
    Ar2 <- Ar2 + 0.5 * pctRec2[j] * pctDv2[j] + Abox2
    if (debug){
      cat("i", i, " j", j, " Ar2", Ar2, " Abox2", Abox2, " pctRec2", pctRec2[j], " pctDv2", pctDv2[j], " class1avg1", class1avg1[j], " dvRunningSum2", dvRunningSum2, "\n")
    }
    # update the running sum
    dvRunningSum2 <- dvRunningSum2 + pctDv2[j]
  }
  # the area for a perfect predictor
  Aperfect <- 1 - 0.5 * sum(tab2[, 2]) / sum(rowSumsTab2)
  if (debug){
    cat("Aperfect",Aperfect,"\n")
  }
  predictivePower <- max(0, (Ar2-0.5) / (Aperfect-0.5) )
  return(predictivePower)
}


PredictivePowerCv <- function(iv, dv, warn.levels=30, debug=FALSE, folds=10, ...){
  if (length(folds)==1) {
    # the number of folds have been specified, set max at 10
    stopifnot(folds>1 & folds<11)
    stopifnot(class(folds)=="numeric")
    # create an index for cross validation folds
    idx <- sample(x=1:folds, size=length(iv), replace=TRUE)
    num.unique <- folds
  } else {
    # require that the input is two or more unique values, up to 11, with same length as iv
    stopifnot(length(folds) == length(iv))
    num.unique <- length(unique(folds))
    stopifnot(num.unique>1 & num.unique<11)
    idx <- folds
  }
  
  # if the independent variable is numeric then convert it to a factor with BinaryCut
  if (class(iv)=="numeric"){
    fiv <- BinaryCut(iv, dv, ...)
  } else if (class(iv)=="factor") {
    fiv <- iv
  } else {
    stop("Class of iv must be numeric or factor, ", class(iv), " provided.")
  }
  
  # loop for each fold of cross validation
  unique.values <- unique(idx)
  predictive.power <- rep(0, num.unique)
  for (i in 1:num.unique){
    # create an index of booleans, 
    idx.fold <- unique.values[i] != idx
    # calculate predictive power
    predictive.power[i] <- PredictivePower(fiv, dv, warn.levels=warn.levels, debug=debug, cv=idx.fold)
  }
  
  return(list(
    predictive.power = predictive.power,
    mean = mean(predictive.power),
    sd = sd(predictive.power),
    # robustness computed with the coefficient of variation (sd over mean)
    robustness = max(0,1-sd(predictive.power)/mean(predictive.power)) ))
}


PredictivePower.numeric <- function(iv, dv, warn.levels=30, cv=NULL, debug=FALSE, ...){
  fiv <- BinaryCut(iv, dv, ...)
  return(PredictivePower.factor(fiv, dv, warn.levels=warn.levels, debug=debug, cv=cv))
}