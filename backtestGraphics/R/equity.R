#' Equity Data from 2005 to 20014.
#' 
#' This data frame contains the information of 50 equities from 2005-05-02 to
#' 2014-04-30. These 50 stocks are randomly chosen from a backtest results of
#' the whole stock markets. 
#' 
#' @format A data frame with 5 variables 
#' \itemize{ 
#'   \item name = The name of this equity.
#'   \item date = trading date.
#'   \item sector = The sector this equity belongs to. There are nine sectors
#'         in total.
#'   \item nmv = The net market value of the equity held on that day.
#'   \item pnl = The adjusted P&L of the equity on that day. P&L is adjusted
#'         according to the bid offer costs.
#' }
#' @docType data
#' @name equity
#' @keywords equity backtest data
NULL