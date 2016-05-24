#' Calculates compounded interest rate
#'
#' Based on loan term, instalment, and the loan amount, this function calculates 
#' the associated compound interest rate.  This function is designed to be 
#' equivalent to the Excel function RATE.  It calculates a fixed interest rate.
#'
#' @param nper Number of periods
#' @param pmt Instalment per period (should be negative)
#' @param pv Present value i.e. loan advance (should be positive)
#' @param fv      Future value i.e. redemption amount
#' @return rate The corresponding compound interest rate required to arrive at an FV of 0
#'
#' 
#' @keywords financial pv pmt rate
#' @seealso \code{\link{PMT}}  \code{\link{PV}} 
#' @family finance
#' @export
#' 
#' @examples
#' RATE(12,-500,3000) # 0.126947 Taken from excel
#' 
#' df<-data.frame(nper=c(12,12),pmt=c(-500,-400),pv=c(3000,3000))
#' RATE(df$nper,df$pmt,df$pv)  # c(0.126947,0.080927) Taken from excel

RATE <- function(nper, pmt, pv, fv = 0) {
    stopifnot(is.numeric(pv), is.numeric(nper), is.numeric(pmt), is.numeric(fv), pv > 0, nper >= 1, pmt < 0)
    rate <- function(nper, pmt, pv, fv) {
        rate1 <- 0.01
        rate2 <- 0.005
        for (i in 1:10) {
            
            pv1 <- PV(rate1, nper, pmt, fv) - pv
            pv2 <- PV(rate2, nper, pmt, fv) - pv
            
            if (pv1 != pv2) {
                newrate <- (pv1 * rate2 - pv2 * rate1)/(pv1 - pv2)
            }
            
            if (abs(pv1) > abs(pv2)) {
                rate1 <- newrate
            } else {
                rate2 <- newrate
            }
        }
        rate1
    }
    rate <- Vectorize(rate)
    return(rate(nper, pmt, pv, fv))
} 
