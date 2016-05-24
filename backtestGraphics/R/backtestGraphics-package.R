#' A package to visualize backtest results
#' 
#' \tabular{ll}{ Package \tab backtestGraphics \cr Type: \tab Package\cr 
#' Version: \tab 0.1-1\cr Date: \tab 2015-07-14\cr License: \tab GPL-3\cr }
#' 
#' The backtestGraphics package creates an interactive graphics interface for 
#' commodity futures portfolios, equity portfolios and credit default swap. The 
#' interface contains three drop-down menus that allow the user to look at 
#' different portfolios, different investment strategies and different sectors 
#' or commodities. Summary statistics of the start date, end date, allocated 
#' capital, average GMV, number of instruments, cumulative P&L, annualized P&L, 
#' annualized volatility, sharpe ratio, best month and worst month are displayed
#' in the summary screen. More details of the top three drawdowns, three best 
#' performers and three worst performers are displayed in another detail screen.
#' An interactive plot for the cumulative P&L or the daily P&L is shown at the 
#' top and another interactive graph for nmv, gmv and number of contracts is 
#' displayed at the bottom.
#' 
#' The package contains a main function \code{backtestGraphics} that takes in a 
#' data frame of backtest results, and returns summary statistics about the data 
#' frame and plots the historical traces of market values and profits. The package
#' also contains many helper functions that perform the calculations and 
#' plotting for the \code{backtestGraphics} function. These helper functions can
#' only be called within the \code{backtestGraphics} function.
#' 
#' The package contains three data frames, \code{commodity}, \code{equity} and 
#' \code{credit}. The \code{commodity} data frame contains the backtest results 
#' for 28 commodities in the futures market. The \code{equity} data frame 
#' contains the backtest results for random 50 stocks. The \code{credit} data 
#' frame is the backtest result for credit default swap. The user can simply 
#' call \code{backtestGraphics} with these data frames. An example will be
#' \code{backtestGraphics(x = commodity)}.
#' 
#' In order to use the user's own data frame, sometimes the user might need to
#' specify the column names of her data frame and pass them into the function.
#' Type \code{?backtestGraphics} for more details.
#' 
#' @name backtestGraphics-package
#' @docType package
NULL