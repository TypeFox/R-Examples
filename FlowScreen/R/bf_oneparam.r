#' One parameter recursive digital filter
#' 
#' This function estimates baseflow.
#' @param discharge Numeric vector of daily flow data
#' @param k Numeric value for the recession constant (dimensionless).
#' @return Returns a numeric vector of the estimated baseflow.
#' @references Eckhardt, K. 2005. How to construct recursive digital filters
#'   for baseflow separation methods. Journal of Hydrology 352: 168-173.
#' @author Paul H. Whitfield
#' @export
#' @examples
#' data(cania.sub.ts)
#' res <- bf_oneparam(cania.sub.ts$Flow, k=0.9)
#' plot(cania.sub.ts$Date, cania.sub.ts$Flow, xlab="", ylab="Q (m3/s)", type="l")
#' points(cania.sub.ts$Date, res, type="l", col="blue")

bf_oneparam <- function(discharge, k){
    bf <- rep(discharge[1],length(discharge))
    for(i in 2:length(discharge)) {
        bf[i] <- (k*bf[i-1]/(2-k)) + ((1-k)*discharge[i]/(2-k))
        if(bf[i] > discharge[i]) bf[i] <- discharge[i]
    }
    return(bf)
}