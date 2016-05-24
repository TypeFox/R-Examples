#############################################################################
##
##  Create guide table.
##  (R version)
##
#############################################################################

make.guidetable.R <- function(params) {

  ## Get number of interval.
  n.ivs <- ncol(params)-1

  ## Allocate memory for guide table.
  gt <- integer(n.ivs)

  ## Compute cumulative sums.
  Acum <- cumsum(params["A.ht", 1:n.ivs])
  A.ht.tot <- Acum[n.ivs]

  ## Compute average area for each interval.
  Astep <- A.ht.tot / n.ivs

  ## Create hash table (guide table).
  sum <- 0
  i <- 1
  for (j in 1:n.ivs) {
    while (is.TRUE(Acum[i] < sum))
      i <- i+1
    if (is.TRUE(i > n.ivs))
      break;
    gt[j] <- as.integer(i-1)   ## This is required for the C version.
    sum <- sum + Astep
  }
  
  ## If there has been an round off error,
  ## we have to complete the guide table.
  if (j<n.ivs) {
    for(jj in j:n.ivs)
      gt[jj] <- n.ivs
  }

  ## Return cumulative areas and hash table.
  return(list(Acum=Acum, gt=gt))
}

## --------------------------------------------------------------------------

## Remark:
## In C arrays are indexed using 0, 1, ..., length(array)-1
## In R arrays are indexed using 1, 2, ..., length(array)
## The guide table is primarily used by the C version.

## --------------------------------------------------------------------------
