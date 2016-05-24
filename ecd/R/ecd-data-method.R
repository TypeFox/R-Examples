#' Read sample data
#' 
#' Read sample data into xts by specifying the symbol.
#' The xts object has two rows: the prices indexed by the dates.
#'
#' @param symbol Character for the symbol of the time series. Default: dji
#'
#' @return The xts object for the time series
#'
#' @keywords timeseries xts sample-data
#'
#' @export
#'
#' @examples
#' dji <- ecd.data()
#' wti <- ecd.data("wti")
### <======================================================================>
"ecd.data" <- function(symbol = "dji")
{
    if(grepl("option", symbol)) {
        stop(paste("Option is not supported for symbol", symbol))
    }
    
    # basic settings
    
    dt <- "Date"
    col_out <- "Close"
    
    c <- .ecd.data_config(symbol)
    if(nrow(c) != 1){
        stop(paste("Unknown symbol", symbol, "for sample data!"))
    }
    
    df <- ecd.read_csv_by_symbol(symbol)
    
    ts <- ecd.df2ts(df, date_format = c$date_format,
                    dt = dt, col_in = c$col_in, col_out = col_out)
    xtsAttributes(ts) <- list(symbol = symbol)
    ts
}
### <---------------------------------------------------------------------->
