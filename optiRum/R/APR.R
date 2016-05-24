#' Calculates the compound interest rate for a loan
#'
#' Based on period interest rate, number of periods, and loan amount, this function calculates 
#' the compound annual interest rate of the loan based on the monthly repayment.
#' It calculates based on a fixed interest rate, FV=0, and charging is 
#' at the end of the period.
#'
#' @param nper    Number of periods - monthly
#' @param pmt     Instalment per period (should be negative)
#' @param pv      Present value i.e. loan advance (should be positive)
#' @param fv      Future value i.e. redemption amount
#' 
#' @return rate   The effective interest rate per year
#'
#' @keywords financial pv pmt apr
#' @seealso \code{\link{RATE}}
#' @family finance
#' @export
#' 
#' @examples
#' # single set of values
#' APR(12,-10,110) 
#' 
#' # vector of values
#' df<-data.frame(nper=c(12,24),pmt=c(-10,-10),pv=c(110,220))
#' APR(df$nper,df$pmt,df$pv)
#' 

APR <- function(nper, pmt, pv, fv = 0) {
    stopifnot(nper >= 1, pmt < 0, pv > 0)
    rate <- ((1 + RATE(nper, pmt, pv, fv))^12) - 1
    return(rate)
} 
