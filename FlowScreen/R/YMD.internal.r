#' Add calendar year, month, and day of year columns
#' 
#' @param TS Output from \code{\link{create.ts}} function.
#' @return Returns a data.frame with year, month, and doy columns appended.
#' @author Jennifer Dierauer

YMD.internal <- function(TS) {
    TS$year <- as.numeric(format(TS$Date, "%Y"))
    TS$month <- as.numeric(format(TS$Date, "%m"))
    TS$doy <- as.numeric(format(TS$Date, "%j"))
    
    return(TS)
}