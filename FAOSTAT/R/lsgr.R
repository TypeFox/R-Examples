##' Least squares growth rate
##'
##' Function for generating the n-period rolling least squares growth
##' rate.
##'
##' Missing values are ommited in the regression. (Will need to check
##' this.)
##'
##' TODO (Michael): There is still some error associated with this
##' function, will need to investigate further. Will need a rule for
##' this, when the fluctuation is large and data are sufficient then
##' take the lsgr, otherwise the geogr.
##'
##' In order to ensure the growth rate calculated is reliable, the
##' following rule are applied.
##' \enumerate{
##'   \item 50\% of the data must be present.
##'   \item The length of the time series must be greater than n.
##' }
##' Otherwise the growth will not be computed.
##'
##' @param x The time series for the growth rate to be calculated
##' @param n The period for the growth to be calculated over.
##' @return The n-period least squares growth rate of the time series
##' @export
##' @examples
##' test.ts = abs(rnorm(100))
##' lsgr(test.ts, 1)
##' lsgr(test.ts, 3)
##' lsgr(test.ts, 10)
##'

lsgr = function(x, n = 1){
    T = length(x)
    if(all(is.na(x))){
        lsgr = rep(NA, T)
        warning("All values are NA")
    } else {
      firstObs = min(which(!is.na(x)))
        if(length(na.omit(x)) < 5){
            stop("Insufficient data for least squares growth rate, use other methods")
        } else {
            if(sum(is.na(x[firstObs:T])) > 0.5 * (T - firstObs + 1)){
                t = 1:(n + 1)
                lsgr = double(T)
                lsgr[1:(firstObs + n - 1)] = NA
                for(i in (firstObs + n):T){
                    tmp = try((exp(coef(rlm(log(x[(i - n):(i)]) ~ t))[2]) - 1) *
                      100, silent = TRUE)
                    if(!inherits(tmp, "try-error")){
                        lsgr[i] = tmp
                    } else {
                        lsgr[i] = NA
                    }
                }
                warning("Over 50% of the data are missing, robust regression is used")
            } else {
                t = 1:(n + 1)
                lsgr = double(T)
                lsgr[1:(firstObs + n - 1)] = NA
                for(i in (firstObs + n):T){
                    tmp = try((exp(coef(lm(log(x[(i - n):(i)]) ~ t))[2]) - 1) * 100, 
                              silent = TRUE)
                    if(!inherits(tmp, "try-error")){
                        lsgr[i] = tmp
                    } else {
                        lsgr[i] = NA
                    }
                }
            }
        }
    }
    lsgr
}
