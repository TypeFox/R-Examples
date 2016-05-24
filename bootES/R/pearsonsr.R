## Daniel Gerlanc and Kris Kirby (2010-2012)
## Function for calculating Pearson's r

calcPearsonsR <- function(vals, freq, grps, lambdas) {
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
  ## @pre: 'lambdas' is not optional and the sum of the lambdas must be 0
  ##  
  ## Details:
  ##   This function is meant to be passed as the 'statistic' argument
  ##   to the 'boot' function. 'freq' should be a frequency vector of
  ##   the type returned by 'boot'
  
  ## Get the integer indices of the different groups.
  grp.idx = split(seq_along(vals), grps, drop=TRUE)   
  grp.nms = names(grp.idx)
  
  ## Calculate the mean and sum-of-squares for the bootstrap samples
  ## from each group.
  ns = means = ss = numeric()
  for (nm in grp.nms) {
    this.grp.idx = grp.idx[[nm]]
    sample.vals  = rep(vals[this.grp.idx], times=freq[this.grp.idx])

    means[nm] = mean(sample.vals) # stats for each group
    ss[nm]    = sum((sample.vals - means[nm])^2)
    ns[nm]    = length(sample.vals)
  }

  ## Put means and lambdas in the same order
  lambdas = lambdas[match(names(lambdas), grp.nms)]
  C       = sum(lambdas * means) # contrast
  sumWts  = sum(lambdas^2 / ns) # sum of the weights
  SSw     = sum(ss)
  
  res =  C / sqrt(C^2 + SSw * sumWts)
  res
}



