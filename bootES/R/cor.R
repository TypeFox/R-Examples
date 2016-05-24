## Daniel Gerlanc and Kris Kirby (2010-2012)
## Function for calculating bootstrap correlations

corBoot <- function(df, i) {
  ## Compute the Pearson Product-Moment correlation r between two variables
  ## using an index vector 'i' to indicate which values should be sampled
  ## 
  ## Args:
  ##  df : A data frame containing 2 numeric columns
  ##  i  : An index vector of the type returned by 'boot' when stype='i'
  
  cor(df[i, 1], df[i, 2])
}
  
calcBootCor <- function(data, freq) {
  ## Compute a bootstrap correlation between two variables
  ## 
  ## Args:
  ##   data: a data frame containing columns 'x' and 'y'
  ##   freq: a frequency vector of length equal to the number of records in 
  ##    'data' indicating how many times an observation should be drawn
  ##     from each sample.
  ##  
  ## Details:
  ##   This function is meant to be passed as the 'statistic' argument
  ##   to the 'boot' function. 'freq' should be a frequency vector of
  ##   the type returned by 'boot' when stype='f'

  sample.idx  = rep(seq_len(nrow(data)), times=freq)
  sample.vals = data[sample.idx, c(1, 2)]
  res         = cor(sample.vals[[1]], sample.vals[[2]])
  res
}

calcBootCorDiff <- function(data, freq, grps) {
  ## Compute the difference in correlation between 2 groups
  ## 
  ## Args:
  ##   data: a data frame containing columns 'x' and 'y'
  ##   freq: a frequency vector of length equal to the number of records in 
  ##    'data' indicating how many times an observation should be drawn
  ##     from each sample.
  ##   grps: a grouping vector of the same length as 'vals'
  ##  
  ## Details:
  ##   This function is meant to be passed as the 'statistic' argument
  ##   to the 'boot' function. 'freq' should be a frequency vector of
  ##   the type returned by 'boot' when stype='f'
  
  ## Get the row indices of the different groups.
  grp.idx = split(seq_len(nrow(data)), grps, drop=TRUE)   

  if (length(grp.idx) != 2)
    stop("There must be only 2 groups!")
  
  ## Calculate the bootstrap correlation for each group.
  cors = numeric()
  for (i in seq_along(grp.idx)) {
    this.grp.idx = grp.idx[[i]]
    sample.idx   = rep(this.grp.idx, times=freq[this.grp.idx])
    sample.vals  = data[sample.idx, c(1, 2)]
    cors[i]      = cor(sample.vals[[1]], sample.vals[[2]])
  }

  res = cors[1] - cors[2]
  res
}

