#' Eckhardt two parameter recursive digital filter 
#' 
#' This function takes vector of discharge data and estimates the baseflow
#' @param discharge vector of daily discharge observations
#' @param a Numeric value.
#' @param BFI Numeric value.
#' @return Returns
#' @author Paul Whitfield
#' @references Eckhardt, K. 2012. Technical note: Analytical sensitivity analysis 
#'   of two parameter recursive digital baseflow separation filter. Hydrology and
#'   Earth System Sciences 16: 451-455.
#' @export
#' @examples
#' data(cania.sub.ts)
#' bf <- bf_eckhardt(cania.sub.ts$Flow, 0.97, 0.8)
#' plot(cania.sub.ts$Date, cania.sub.ts$Flow, type="l")
#' points(cania.sub.ts$Date, bf, type="l", col="blue")

bf_eckhardt <- function(discharge, a, BFI){
    bf <- rep(discharge[1],length(discharge))
    for(i in 2:length(discharge)) {
        bf[i] <-(((1 - BFI)* a* bf[i-1]) + ((1-a)* BFI* discharge [i])) /(1- a*BFI)
        if(bf[i] > discharge[i]) bf[i] <- discharge[i]
    }
    return(bf)
}