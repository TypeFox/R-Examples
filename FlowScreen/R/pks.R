#' Get the flow peaks over a threshold
#' 
#' This function finds the flow peaks over a user defined threshold and declusters to remove
#' dependent peaks.
#' @param TS output from \code{\link{create.ts}} containing a data.frame of flow
#'   time series
#' @param Dur numeric value of the minimum number of days between peaks
#' @param Qmax numeric value for peaks over threshold quantile.  
#'   Default is 0.95.
#' @details Peaks Over Threshold values are calcuated as mean daily streamflow (m3/s)
#'   minus the threshold streamflow value (m3/s) defined by the 
#'   input quantile (Qmax). Peaks are identified with \code{\link[evir]{pot}} and 
#'   the minimum inter-event duration (Dur) is applied by 
#'   \code{\link[evir]{decluster}}.
#' @return Returns a numeric vector of peaks of threshold values in m3/s. The "times"
#'   attribute contains the date by calendar year, and the "names" attribute contains
#'   the hydrologic year and hydrologic day of year, e.g., 2012 55.
#' @author Jennifer Dierauer
#' @export
#' @examples
#' data(cania.sub.ts)
#' res <- pks(cania.sub.ts)
#' res2 <- screen.metric(res, "Peak Over Threshold (m3/s)")

pks <- function(TS, Dur=5, Qmax=0.95) {
    
    MyThreshold <- stats::quantile(TS$Flow, Qmax, na.rm=TRUE)
    
    Flow<-TS$Flow
    attr(Flow, "times") <- TS$Date
    attr(Flow, "names") <- paste(TS$hyear, TS$hdoy, sep=" ")
    
    ##Find Peaks Over Threshold
    Peaks <- evir::pot(Flow, MyThreshold, picture=F)
    
    ##Remove minor peaks using input "duration" between peaks
    Peaks <- evir::decluster(Peaks$data, Dur, picture=F)
    Peaks <- Peaks - MyThreshold
    return(Peaks)
}
