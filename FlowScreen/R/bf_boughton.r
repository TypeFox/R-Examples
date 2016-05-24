#' Boughton recursive digital filter
#' 
#' This function estimates baseflow
#' @param discharge Nnumeric vector of daily flow data
#' @param k Numeric value of the recession constant (dimensionless).
#' @param C Numeric value of the partitioning factor (dimensionless).
#' @return Returns a numeric vector of the estimated baseflow.
#' @references Boughton, WC. 1993. A hydrograph-based model for estimating the 
#'   water yield of ungauged catchments.In Hydrology and Water Resources Symposium,
#'   Institution of Engineers Australia, Newcastle, NSW; 317-324.
#' @author Paul H. Whitfield
#' @export
#' @examples
#' data(cania.sub.ts)
#' res <- bf_boughton(cania.sub.ts$Flow, k=0.9, C=0.1)
#' plot(cania.sub.ts$Date, cania.sub.ts$Flow, xlab="", ylab="Q (m3/s)", type="l")
#' points(cania.sub.ts$Date, res, type="l", col="blue")

bf_boughton <- function(discharge, k, C){
    bf <- rep(discharge[1],length(discharge))
    for(i in 2:length(discharge)) {
        bf[i] <- (k*bf[i-1]/(1+C)) + (C*discharge[i]/(1+C))
        if(bf[i] > discharge[i]) bf[i] <- discharge[i]
    }
    return(bf)
}