#' Read option data csv
#'
#' Read option data csv into dataframe, enriched with
#' \code{Date, expiration_date, days}.
#'
#' @param symbol character, option data symbol
#'
#' @return dataframe
#'
#' @keywords data
#'
#' @author Stephen H-T. Lihn
#'
#' @export
#'
#' @examples
#'
#' df <- ecop.read_csv_by_symbol("spxoption2")
#'
### <======================================================================>
"ecop.read_csv_by_symbol" <- function(symbol)
{
    df <- ecd.read_csv_by_symbol(symbol)
    date_format = "%Y%m%d"
    df$Date <- as.Date(as.character(df[,"TRADE_DT"]), date_format)
    df$expiration_date <- as.Date(as.character(df[,"EXPR_DT"]), date_format)
    df$days <- as.numeric(difftime(df$expiration_date, df$Date, units="days"))
    df
}
### <---------------------------------------------------------------------->
