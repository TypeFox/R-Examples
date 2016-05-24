#' Daily S&P 500 data from Jan 3, 1950 to present
#'
#' @description
#' The "usual" S&P 500 index, without dividends reinvested.
#'
#' @param tr (bool) choose Total Return (tr=TRUE) or not (default tr=FALSE)
#' Total Return data is available only from 1988; non-dividend-reinvested data is available from 1950.
#'
#' @return A data.frame with the daily data
#'
#' @references
#' Yahoo Finance
#'
#' @details The columns of the data.frame returned
#' \itemize{
#'     \item \bold{Date}
#'     \item \bold{Open}
#'     \item \bold{High}
#'     \item \bold{Low}
#'     \item \bold{Close}
#'     \item \bold{Volume}
#'     \item \bold{Adj.Close}
#'     }
#'
#' @importFrom utils read.table
#'
#' @author George Fisher
#'
#' @examples
#' sp500_idx <- SP500()
#' head(sp500_idx)
#' tail(sp500_idx)
#'
#' @export
SP500 <- function(tr=FALSE) {
    if (tr) {
        URL = "http://ichart.finance.yahoo.com/table.csv?s=%5ESP500TR"
    } else {
        URL = "http://ichart.finance.yahoo.com/table.csv?s=%5EGSPC"
    }
    table <- read.table(URL,
                        header = TRUE,sep=",",stringsAsFactors = FALSE)
    table$Date <- as.Date(table$Date, "%Y-%m-%d")
    return(table)
}
