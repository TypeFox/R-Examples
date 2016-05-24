#' Calculates the present value
#'
#' Based on period interest rate, number of periods, and instalment, this function calculates 
#' the present value of the loan such that it would be paid off fully at the end of the loan.
#' This function is designed to be equivalent to the Excel function PV. 
#' It calculates based on a fixed interest rate, FV=0 and charging is 
#' at the end of the period. Response is rounded to 2dp
#'
#' @param rate The nominal interest rate per period (should be positive)
#' @param nper Number of periods
#' @param pmt Instalment per period (should be negative)
#' @param fv      Future value i.e. redemption amount
#' @return pv Present value i.e. loan advance (should be positive)
#'
#' @keywords financial pv pmt
#' @seealso \code{\link{PMT}}  \code{\link{RATE}} 
#' @family finance
#' @export
#' 
#' @examples
#' PV(0.1,12,-10) # 68.14 Taken from excel
#' 
#' df<-data.frame(rate=c(.1,.1),nper=c(12,24),pmt=c(-10,-15))
#' PV(df$rate,df$nper,df$pmt)  # c(68.14,134.77) Taken from excel


PV <- function(rate, nper, pmt, fv = 0) {
    stopifnot(is.numeric(rate), is.numeric(nper), is.numeric(pmt), is.numeric(fv), rate > 0, rate < 1, nper >= 1, pmt < 0)
    
    pvofregcash <- -pmt/rate * (1 - 1/(1 + rate)^nper)
    pvoffv <- fv/((1 + rate)^nper)
    
    return(round(pvofregcash - pvoffv, 2))
} 
