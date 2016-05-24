################################################################################
##
## $Id: bt.spread.R 1300 2008-08-27 21:01:11Z zhao $
##
## Calculates the spreads and confidence intervals for a given data frame
##
################################################################################

## NOTE: the confidence interval is a hack. We assume that the spread
## is a sort of weighted mean calculation in which the weights are 1
## for the long quantile and 1 for the short quantile.

## "m" is a 2-dimensional array of means for each in.var or
## in.var/by.var.  Normally accessed through the 5th dimension of the
## "results" slot (object@results[, , , , "means"]).

## "n" is a 2-dimensional array that contains the number of
## observations (count) for each in.var or in.var/by.var combination.
## Normally accessed through the 5th dimension of the "results" slot
## (object@results[, , , , "counts"]).

## "sd" is the standard deviation of a specific measure of return.
## Normally stored in ret.stats[ret.var, "sd"] where ret.stats is a
## slot of "backtest" and ret.var is the return variable for which we
## want the standard deviation

.bt.spread <- function(m, n, sd){

  ## subtracts the mean of highest the quantile from the mean of the
  ## lowest quantile

  spread <- m[ ,dim(m)[2]] - m[ ,1]
  
  ## calculates standard error

  se <- sd[1]/sqrt(n[ ,1] + n[ ,dim(n)[2]])
  
  ## 95.46% of values fall w/in 2 standard deviations

  ci.low <- spread - (2 * se)
  ci.high <- spread + (2 * se)

  result <- cbind(spread, ci.low, ci.high)
  
  result <- array(result, dim = dim(result), dimnames =
                  list(dimnames(m)[[1]], c("spread", "CI(low)", "CI(high)")))
  
  result
  
}
