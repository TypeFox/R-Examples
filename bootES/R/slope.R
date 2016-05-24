## Daniel Gerlanc and Kris Kirby (2010-2012)
## Functions for calculating slope

calcSlope <- function(vals, freq, grps, lambdas) {
  ## Compute Pearson's r for one or more groups
  ## 
  ## Args:
  ##   vals: a numeric vector of values
  ##   freq: a frequency vector of length equal to the number of records in 
  ##    'data' indicating how many times an observation should be drawn
  ##     from each sample.
  ##   grps: a grouping vector of the same length as 'vals'
  ##   lambdas: an optional named vector of contrast weights
  ##  
  ## Details:
  ##   This function is meant to be passed as the 'statistic' argument
  ##   to the 'boot' function. 'freq' should be a frequency vector of
  ##   the type returned by 'boot'
  
  ## Get the integer indices of the different groups.
  grp.idx = split(seq_along(vals), grps, drop=TRUE)   
  grp.nms = names(grp.idx)
  
  ## Calculate the mean for the bootstrap samples from each group.
  means = numeric()
  for (nm in grp.nms) {
    this.grp.idx = grp.idx[[nm]]
    sample.vals  = rep(vals[this.grp.idx], times=freq[this.grp.idx])
    means[nm] = mean(sample.vals) # stats for each group
  }

  ## Put means and lambdas in the same order
  lambdas = lambdas[match(names(lambdas), grp.nms)]
  res     = sum(lambdas * means) # contrast
  res
}


