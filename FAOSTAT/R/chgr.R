##' Absolute change between the year
##'
##' Function for generating the n-period absolute change
##'
##' In order to ensure the change calculated is reliable, the
##' following rule are applied.
##' \enumerate{
##'   \item 50\% of the data must be present.
##'   \item The length of the time series must be greater than n
##' }
##' Otherwise the growth will not be computed.
##'
##' @param x The time series for the change to be calculated.
##' @param n The period for the growth to be calculated over.
##' @return The n-period change of the time series.
##' @export
##' @examples
##' test.ts = abs(rnorm(100))
##' chgr(test.ts, 1)
##' chgr(test.ts, 3)
##' chgr(test.ts, 10)

chgr = function(x, n = 1){
  T = length(x)
  if(sum(is.na(x)) == T){
    chgr = rep(NA, T)
    warning("All values are NA")
  } else {
    firstObs = ifelse(any(is.na(x)), min(c(1, which(!is.na(x)))), 1)
    if(n > T - firstObs - 1){
      chgr = rep(NA, T)
      warning("Time series does not have sufficient values")
    } else {
      if(sum(is.na(x[firstObs:T])) > 0.5 * (T - firstObs + 1)){
        chgr = rep(NA, T)
        warning("Over 50% of the data are missing")
      } else {
        chgr = double(T)
        chgr[1:(firstObs + n - 1)] = NA
        chgr[(firstObs + n):T] = diff(x[(firstObs):T], n)
      }
    }
  }
  chgr
}
