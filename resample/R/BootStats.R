# Copyright 2014 Google Inc. All rights reserved.
#
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file or at
# http://opensource.org/licenses/BSD-3-Clause

.BootStats <- function(x) {
  # return a data frame with basic summary statistics
  if(is.null(x$replicates))
    return(NULL)
  meanReplicates <- colMeans(x$replicates)
  data.frame(row.names = names(x$observed),
             Observed = x$observed,
             SE = colStdevs(x$replicates),
             Mean = meanReplicates,
             Bias = meanReplicates - x$observed)
}

.PermutationStats <- function(x, alternative,
                              tolerance = .Machine$double.eps ^ 0.5) {
  # return a data frame with basic summary statistics
  #
  # Args:
  #   x: a permutationTest or permutationTest2 object
  #   alternative: character vector: "greater", "less", or "two.sided"
  #      This is replicated to length = number of statistics
  #   tolerance: differences smaller than tolerance (on either absolute
  #      or relative scales) are considered equal

  if(is.null(x$replicates))
    return(NULL)
  if(any(!is.element(alternative, c("greater", "less", "two.sided"))))
    warning("alternative is not a recognized value")
  alternative <- rep(alternative, length = x$p)

  R <- x$R
  o <- x$observed
  r <- x$replicates
  tol <- tolerance * max(1, colStdevs(r))
  pValueLess <- (1 + colSums(r <= rep(o + tol, each = R))) / (R+1)
  pValuePlus <- (1 + colSums(r >= rep(o - tol, each = R))) / (R+1)
  pValue <- rep(NA, x$p)
  pValue[alternative == "less"] <- pValueLess[alternative == "less"]
  pValue[alternative == "greater"] <- pValuePlus[alternative == "greater"]
  pValue[alternative == "two.sided"] <- pmin(1,
           2 * pmin(pValueLess, pValuePlus))[alternative == "two.sided"]

  data.frame(row.names = names(o),
             Observed = o,
             Mean = colMeans(x$replicates),
             Alternative = alternative,
             PValue = pValue)
}

.JackknifeStats <- function(x) {
  # return a data frame with basic summary statistics
  meanReplicates <- colMeans(x$replicates)
  n <- x$n
  data.frame(row.names = names(x$observed),
             Observed = x$observed,
             SE = (n-1)/sqrt(n) * colStdevs(x$replicates),
             Mean = meanReplicates,
             Bias = (n-1) * (meanReplicates - x$observed))
}
