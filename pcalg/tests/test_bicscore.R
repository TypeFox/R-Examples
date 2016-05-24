#' Tests the calculation of BIC and MLE
#' 
#' @author Alain Hauser
#' $Id: test_bicscore.R 224 2014-02-06 14:23:53Z alhauser $

cat("Testing the calculation of the BIC score... ")

library(pcalg)

load("test_bicscore.rda") # in directory tests/ i.e., typically *not* installed

## Tolerance for numerical comparison
tol <- sqrt(.Machine$double.eps)

for (cpp in c(FALSE, TRUE)) {
  score <- new("GaussL0penIntScore", 
      targets = gauss.targets, 
      target.index = gauss.target.index, 
      data = gauss.data,
      use.cpp = cpp)
  
  # print(score$pp.dat)

  if (any(score$pp.dat$data.count != 1000))
    stop("The number of non-interventions are not calculated correctly.")

  if (any(score$pp.dat$scatter.index != 1:5))
    stop("The indices of the scatter matrices are not calculated correctly.")

  for (i in 1:5)
    if (!isTRUE(all.equal(score$pp.dat$scatter[[score$pp.dat$scatter.index[i]]][1:5, 1:5], 
        gauss.scatter[[i]], 
        tolerance = tol)))
      stop("The scatter matrices are not calculated correctly.")

  for (i in 1:5)
    if (!isTRUE(all.equal(gauss.loc.score[[i]], 
        score$local.score(i, gauss.parents[[i]]), 
        tolerance = tol)))
      stop("The local score is not calculated correctly.")
  
  # print(lapply(1:5, function(i) score$local.mle(i, gauss.parents[[i]], DEBUG.LEVEL = 3)))
  
  for (i in 1:5) {
    local.mle <- score$local.mle(i, gauss.parents[[i]])
    if (length(local.mle) != length(gauss.mle[[i]]) ||
        !isTRUE(all.equal(gauss.mle[[i]], local.mle, 
        tolerance = tol)))
      stop("The local MLE is not calculated correctly.")
  }
}

cat("Done.")
