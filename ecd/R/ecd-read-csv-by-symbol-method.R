#' Read csv file of sample data
#' 
#' This is a helper utility to read sample csv file into data frame.
#' The main use for external users is to read the option data
#' since it has a different format than other price timeseries data.
#'
#' @param symbol Character for the symbol of the time series. Default: dji
#'
#' @return The data.frame object
#'
#' @keywords timeseries sample-data
#'
#' @export
#'
#' @importFrom utils read.csv
#'
#'
#' @examples
#' dji <- ecd.read_csv_by_symbol("dji")
#' spx <- ecd.read_csv_by_symbol("spxoption2")
#'
### <======================================================================>
"ecd.read_csv_by_symbol" <- function(symbol = "dji")
{
    c <- .ecd.data_config(symbol)
    if(nrow(c) != 1){
        stop(paste("Unknown symbol", symbol, "!\n"))
    }
    file <- paste(c$symbol, "_archive_", c$cols, ".csv", sep="")
    .ecd.read_file(file)
}
### <---------------------------------------------------------------------->
".ecd.read_file" <- function(filename)
{
    # find the sample data location inside package
    locate <- function (filename) {
        f1 <- system.file("extdata", filename, package = "ecd")
        if(length(f1) > 0 & file.exists(f1)) return(f1)
        # during development, this is where it is!
        f2 <- system.file("inst", "extdata", filename, package = "ecd")
        if(length(f2) > 0 & file.exists(f2)) return(f2)
        return("")
    }
    
    # regular file
    f <- locate(filename)
    if (length(f) > 0 & file.exists(f)) {
        return(read.csv(f, header=TRUE))    
    }
    
    # try zip file (arranged to save storage)
    f <- locate(paste(filename, "zip", sep="."))
    if (length(f) > 0 & file.exists(f)) {
        df <- read.csv(unz(f, filename), header=TRUE)
        return(df)
    }
    stop(paste("Failed to locate file for", filename))
}
### <---------------------------------------------------------------------->

