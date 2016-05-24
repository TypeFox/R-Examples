#' Utility to standardize timeseries from data.frame to xts
#' 
#' This utility converts the df input to an xts object with columns
#' and statistics required for the fitting/plot utility in the ecd package.
#' The require columns are Date, Close, logr. This utility can also be used
#' to convert the input from Quandl.
#'
#' @param df          Data.frame of the time serie
#' @param date_format Character, date format of the input date column. 
#'                    It can be NULL to indicate no date conversion is needed.
#'                    Default: "\code{\%m/\%d/\%Y}".
#' @param dt          Character, the name of the input date column. Default: "Date"
#' @param col_in      Character, the name of the input closing price column. Default: "Close"
#' @param col_out     Character, the name of the output closing price column. Default: "Close"
#' @param do.logr     logical, if \code{TRUE} (default), produce xts object of logr; otherwise, just the \code{col_out} column.
#' @param rnd.zero    numeric, a small random factor to avoid an unreal peak of zero log-returns.
#'
#' @return The xts object for the time series
#'
#' @keywords timeseries xts sample-data
#'
#' @export
#'
#' @importFrom stats rnorm
#'
#' @examples
#' \dontrun{
#' ecd.df2ts(df)
#' }
### <======================================================================>
"ecd.df2ts" <- function (df, date_format = "%m/%d/%Y", 
                         dt = "Date",
                         col_in = "Close",
                         col_out = "Close",
                         do.logr = TRUE,
                         rnd.zero = 0.01)
{
    dates <- if (is.null(date_format)) as.Date(df[,dt]) else as.Date(df[,dt], date_format)
    prices <- as.numeric(df[,col_in])
    ts <- xts(prices, dates)
    colnames(ts) <- c(col_out)
    if (!do.logr) return(ts)
    
    # derive log returns
    ts$logr <- diff(log(ts))
    ts <- ts[!is.na(ts$logr)] # remove NA
    sd <- rnorm(length(ts$logr), 0, sd(ts$logr))
    ts$logr <- ts$logr+ifelse(ts$logr != 0, 0, sd)
    ts
}
### <---------------------------------------------------------------------->
