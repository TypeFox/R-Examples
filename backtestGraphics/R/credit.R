#' Credit Default Swap Data from 2007 to 2009.
#' 
#' This data frame contains the information of 260 credit default swaps (CDS) from
#' 2007-01-02 to 2009-12-31. The data sets is actually a combination of CDS backtest
#' conducted under daily, weekly, monthly, and quartrly trading frequency. The strategy
#' column contains the trading frequency.
#' 
#' @format A data frame with 7 variables
#' \itemize{
#'    \item name = The name of each credit default swap (CDS).
#'    \item date = The trading date.
#'    \item sector = The sector the CDS belongs to.
#'    \item strategy = The trading strategy of the credit default swap, including
#'          "daily", "weekly", "monthly" and "quarterly".
#'    \item gmv = The gross market value of the CDS held on that day.
#'    \item nmv = The net market value of the CDS held on that day.
#'    \item pnl = The P&L value (adjusted) of the CDS on that day.
#' }
#' @docType data
#' @name credit
#' @keywords CDS backtest data
NULL