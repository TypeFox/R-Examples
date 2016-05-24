## Daniel Gerlanc and Kris Kirby (2010-2012)

## Functions for calculating unstandardized means and differences in
## unstandardized means. 'calcUnstandardizedMean' is the most
## flexible, allowing the means to be weighted and one or more
## distinct samples. 'meanDiffBoot' and 'meanUnweightedBoot' are more
## specialized, and duplicate some of the functionality of
## 'calcUnstandardizedMean'.

calcUnstandardizedMean <- function(vals, freq, grps, lambdas) {
  ## Compute an unstandardized mean of one or more groups
  ## 
  ## Args:
  ##   vals: a numeric vector of values
  ##   freq: a frequency vector of length equal to the number of records in 
  ##    'data' indicating how many times an observation should be drawn
  ##     from each sample.
  ##   grps: an optional grouping vector of the same length as 'vals'
  ##   lambdas: an optional named vector of weights to scale the means
  ##
  ## Details:
  ##   This function is meant to be passed as the 'statistic' argument
  ##   to the 'boot' function. 'freq' should be a frequency vector of
  ##   the type returned by 'boot' If no 'grps' argument is passed, then
  ##   all values are treated as if they are in the same group. If no 'lambdas'
  ##   argument is passed, then each group is assigned the same weight
  ##
  ## Note:
  ##   To calculate a difference in means, set the lambdas to +1, -1
  ##
  ## Returns:
  ##   The bootstrap, unstandardized mean

  ## Handle a missing 'grps' argument by assigning all values to the same group
  if (missing(grps))
    grps = rep.int(1, length(vals))
  
  ## Get the integer indices of the different groups.
  grp.idx = split(seq_along(vals), grps, drop=TRUE)   
  grp.nms = names(grp.idx)
  n.grps  = length(grp.nms)

  ## Default behavior for missing lambdas is to give each group an equal weight
  if (missing(lambdas))
    lambdas = structure(rep.int(1 / n.grps, n.grps), names=grp.nms)
  
  ## Calculate the mean for each group
  means = numeric()
  for (nm in grp.nms) {
    this.grp.idx = grp.idx[[nm]]
    sample.vals  = rep(vals[this.grp.idx], times=freq[this.grp.idx])
    means[nm]    = mean(sample.vals)
  }

  ## Put means and lambdas in the same order
  lambdas = lambdas[names(means)]
  res   = sum(lambdas * means)
}

meanUnweightedBoot <- function(vals, freq, grps) {
  ## Computes the mean of the means of multiple samples
  ## 
  ## Args:
  ##   vals: a numeric vector of values
  ##   freq: a frequency vector of the same length as 'vals'
  ##     indicating how many times an observation should be drawn
  ##     from each sample.
  ##   grps: a required grouping vector of the same length as 'vals'
  ##
  ## Details:
  ##   This function is meant to be passed as the 'statistic' argument
  ##   to the 'boot' function. 'freq' should be a frequency vector of
  ##   the type returned by 'boot' 
  ##
  ## Returns:
  ##   The mean of the bootstrap means of multiple samples
  
  ## Error Handling
  if (missing(grps))
    stop("'grps' argument cannot be missing.")

  ## Get the indices of the values corresponding to the different samples.
  idx = split(seq_along(vals), grps, drop=TRUE)
  n   = length(idx)
  
  means = numeric()
  for (i in idx) {
    sample.vals = rep(vals[i], times=freq[i])
    means       = append(means, mean(sample.vals))
  }
  res = mean(means)
  res
}

meanDiffBoot <- function(vals, freq, grps) {
  ## Computes the difference in means between 2 samples
  ## 
  ## Args:
  ##   vals: a numeric vector of values
  ##   freq: a frequency vector of the same length as 'vals'
  ##     indicating how many times an observation should be drawn
  ##     from each sample.
  ##   grps: a required grouping vector of the same length as 'vals'
  ##
  ## Details:
  ##   This function is meant to be passed as the 'statistic' argument
  ##   to the 'boot' function. 'freq' should be a frequency vector of
  ##   the type returned by 'boot' 
  ##
  ## Returns:
  ##   The difference in means between 2 samples

  ## Error Handling
  if (missing(grps))
    stop("'grps' argument cannot be missing.")

  ## Get the indices of the values corresponding to the different samples.
  idx = split(seq_along(vals), grps, drop=TRUE)
  n   = length(idx)
  if (n != 2)
    stop("There must be 2 groups for an estimation of differences in means. Only ", n,
         " groups were found.")

  means = numeric()
  for (i in idx) {
    sample.vals = rep(vals[i], times=freq[i])
    means       = append(means, mean(sample.vals))
  }
  res <- means[1] - means[2]
  res
}
