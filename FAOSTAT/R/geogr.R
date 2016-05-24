##' Geometric growth rate
##'
##' Function for generating the n-period rolling geometric growth rate.
##'
##' In order to ensure the growth rate calculated is reliable, the
##' following rule are applied.
##' \enumerate{
##'   \item 50\% of the data must be present.
##'   \item The length of the time series must be greater than n
##' }
##' Otherwise the growth will not be computed.
##'
##' @param x The time series for the growth rate to be calculated.
##' @param n The period for the growth to be calculated over.
##' @return The n-period geometric growth rate of the time series.
##' @export
##' @examples
##' test.ts = abs(rnorm(100))
##' geogr(test.ts, 1)
##' geogr(test.ts, 3)
##' geogr(test.ts, 10)
geogr = function(x, n = 1){
  T = length(x)
  if(sum(is.na(x)) == T){
    geogr = rep(NA, T)
    warning("All values are NA")
  } else {
    firstObs = ifelse(any(is.na(x)), min(which(!is.na(x))), 1)
    if(n > T - firstObs - 1){
      geogr = rep(NA, T)
      warning("Time series does not have sufficient values")
    } else {
      geogr = double(T)
      geogr[1:(firstObs + n - 1)] = NA
      geogr[(firstObs + n):T] = ((x[(firstObs + n):T]/
                                  x[firstObs:(T - n)])^(1/n) - 1) * 100 
    }
  }
  geogr
}
