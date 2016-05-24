# Returns the test sets indices for k-Fold cross-validation. Stratified sampling is used.
# y: the response vector, for stratifying
# k: number of folds
# returns: a list of length k, containing the test-set indices.
CreateFolds <- function(y, k=10) {
  n <- length(y)
  if (n == 0)
    stop('response length is zero')

  uniqY <- unique(y)
  if (!is.factor(y) && length(y) / length(uniqY) >= k) {
# Intepret the integer-valued y as class labels. Stratify if the number of class labels is <= 5.
    y <- factor(y)
  } else if (is.numeric(y)) { 
# 5-stratum Stratified sampling
    if (n >= 5 * k) {
      breaks <- unique(quantile(y, probs=seq(0, 1, length.out=5)))
      y <- as.integer(cut(y, breaks, include.lowest=TRUE))
    } else 
      y <- rep(1, length(y))
  }
  
  sampList <- tapply(seq_along(y), y, SimpleFolds, k=k, simplify=FALSE)
  list0 <- list()
  length(list0) <- k
  samp <- Reduce(function(list1, list2) {
                   mapply(c, list1, list2, SIMPLIFY=FALSE)
  }, sampList, list0)

  return(samp)
}

# Simple k-fold test-set samples.
# Input a set of SAMPLES
# Returns: a list of length k, containing the SAMPLES.
SimpleFolds <- function(yy, k=10) {
  if (length(yy) > 1)
    allSamp <- sample(yy)
  else
    allSamp <- yy

  n <- length(yy)
  nEach <- n %/% k
  samp <- list()
  length(samp) <- k
  for (i in seq_along(samp)) {
    if (nEach > 0)
      samp[[i]] <- allSamp[1:nEach + (i - 1) * nEach]
    else
      samp[[i]] <- numeric(0)
  }
  restSamp <- allSamp[seq(nEach * k + 1, length(allSamp), length.out=length(allSamp) - nEach * k)]
  restInd <- sample(k, length(restSamp))
  for (i in seq_along(restInd)) {
    sampInd <- restInd[i]
    samp[[sampInd]] <- c(samp[[sampInd]], restSamp[i])
  }

  return(samp)
}
