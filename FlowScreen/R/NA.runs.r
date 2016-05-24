#' Missing data runs for daily time series.
#' 
#' This function takes a data.frame from create.ts and returns a data.frame of 
#' missing data runs.
#' @param TS output from \code{\link{create.ts}} containing a data.frame of flow
#'   time series
#' @author Jennifer Dierauer
#' @return Returns a data.frame with the following columns:
#'   \itemize{
#'     \item Start - Date of the start of the missing data period
#'     \item End - Date of the end of the missing data period
#'     \item Duration - number of days in the missing data period
#'   }
#' @seealso \code{\link{create.ts}} to create input, \code{\link{NA.sum}} to sum the
#'   the missing data occurrences by year or month.
#' @export
#' @examples
#' data(caniapiscau)
#' cania.sub <- caniapiscau[300:1200,]
#' cania.ts <- create.ts(cania.sub)
#' res <- NA.runs(cania.ts)


NA.runs <- function(TS) {
    total <- sum(is.na(TS$Flow))
    if (("year" %in% colnames(TS))==FALSE) {
        MyData <- create.ts(TS)
    }
    
    temp <- subset(TS, !is.na(TS$Flow))
    nlength <- length(temp$Flow)
    starts <- as.Date(NA)
    ends <- as.Date(NA)
    for (i in 2:nlength) {
        mlength <- temp$Date[i] - temp$Date[i-1]
        if (mlength > 1) {
            starts <- c(starts, as.Date(temp$Date[i-1]))
            ends <- c(ends, as.Date(temp$Date[i]))
        }
    }
    miss <- data.frame(Start=starts[-1], End=ends[-1])
    miss$Duration <- as.numeric(miss$End-miss$Start)
    
    return(miss)
    
}