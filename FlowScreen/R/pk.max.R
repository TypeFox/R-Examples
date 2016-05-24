#' Annual maximum series
#' 
#' This function returns the annual maximum series from a daily streamflow time
#' series.
#' @param TS output from \code{\link{create.ts}} containing a data.frame of 
#'   the daily streamflow time series
#' @return Returns a numeric vector containing the annual maximum flow (m3/s) 
#'   series, by hydrologic year. The "times" attribute contains the hydrologic year 
#'   for each element in the vector.
#' @author Jennifer Dierauer
#' @seealso See \code{\link{create.ts}} to format the input flow series.
#' 
#'   See \code{\link{pk.max.doy}} to find the day of year for each annual
#'   maximum flow event.
#' @export
#' @examples
#' data(cania.sub.ts)
#' res <- pk.max(cania.sub.ts)
#' res2 <- screen.metric(res, "Q (m3/s)")

pk.max <- function(TS) {
    year <- as.factor(TS$hyear)
    maxseries <- tapply(TS$Flow, year, max, na.rm=T)
    attr(maxseries, "times") <- as.numeric(as.character(unique(TS$hyear)))
    return(maxseries)
}