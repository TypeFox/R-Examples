compCDF <- function(data, weights, 
                    x=seq(min(data, na.rm=TRUE), max(data, na.rm=TRUE), len=250), 
                    comp=1:NCOL(weights), makeplot=TRUE, ...) {
  if (NROW(weights) != NROW(data)) {
    stop("data and weights arguments must have same number of rows")
  }
    
  # First, normalize the weights so the sum of each column is 1/NCOL(data)
  weights <- t(t(weights) / (NCOL(data) * colSums(weights)))
  # Next, give a binomial count for each row of the data and for each x
  f <- function(row, cutpt) colSums(outer(row, cutpt, "<="), na.rm = TRUE)
  bc <- apply(data, 1, f, x)
  # bc is a length(x) by n matrix; each column should be multiplied by
  # the appropriate weight(s) and then the rows summed to give the 
  # unnormalized cdf estimates. This is just a matrix product.
  cdfs <- bc %*% weights[,comp,drop=FALSE]
  
  if(makeplot) {
    plot(range(x), 0:1, type="n", ...)
    for (i in 1:length(comp)) {
      lines(x, cdfs[,comp[i]], lty=i, ...)
    }
  }
  t(cdfs)
}

