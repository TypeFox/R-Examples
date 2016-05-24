#' Get the partial duration series for streamflow droughts
#' 
#' This function returns the partial duration series for streamflow droughts 
#' based on a moving window quantile threshold.
#' @param TS output from \code{\link{create.ts}} containing a data.frame of flow
#'   time series
#' @param Qdr Numeric value of the drought threshold quantile.  Default is 0.2.
#' @param WinSize Numeric value specifying the size of the moving window in
#'   days.  Default is 30. 
#' @details This function defines a daily streamflow threshold and finds the 
#'   partial duration series of streamflow droughts.  Drought events are 
#'   identified in the daily streamflow time series with the threshold level approach.
#'   In this function, the threshold is 
#'   defined by a moving quantile, where daily threshold values are based on
#'   the 80th percentile of the flow duration curve from a 30-day moving window 
#'   (Beyene et al. 2014). With this method, every day of the year has a different
#'   threshold based on the streamflow measured on the day, the 15 days before 
#'   the day, and the 15 days after the day. The size of the moving window can be 
#'   modified with the WinSize argument, and the percentile can be modified with the
#'   Qdr argument.
#' @return Returns the input TS data.frame with "Thresh" and "BelowThresh" 
#'   columns appended. The Thresh column contains the daily flow threshold, and 
#'   the BelowThresh column is a binary indicating whether the flow on each day 
#'   was below the drought threshold.
#' @author Jennifer Dierauer
#' @references Beyene, B.S., Van Loon, A.F., Van Lanen, H.A.J., Torfs, P.J.J.F., 2014. 
#'   Investigation of variable threshold level approaches for hydrological drought 
#'   identification. Hydrol. Earth Syst. Sci. Discuss. 11, 12765-12797. 
#'   http://dx.doi.org/10.5194/hessd-11-12765-2014.
#' @seealso See \code{\link{create.ts}} to format the input flow series.
#' 
#'   See \code{\link{mqt}} to return only the daily moving quantile threshold.
#'   
#'   See \code{\link{dr.events}} to pool drought events, remove minor events,
#'   and calculate metrics.
#'   
#'   See \code{\link{dr.seas}} to calculate metrics for streamflow droughts 
#'   that start in a specific month or months.
#' @export
#' @examples
#' data(cania.sub.ts)
#' pds <- dr.pds(cania.sub.ts)
#' pds <- subset(pds, pds$BelowThresh==TRUE)
#' 
#' # plot the flow time series with black and the drought events in red
#' plot(cania.sub.ts$Date, cania.sub.ts$Flow, ylab="m3/s", xlab="", type="l")
#' points(pds$Date, pds$Flow, pch=19, cex=0.7, col="red")


## input Flow and Date vectors, plus threshold value, e.g. 0.20 after Van Loon and Laaha 2015


dr.pds <- function(TS, Qdr=0.2, WinSize=30) {
    
    doy <- as.factor(TS$doy)
    temp <- TS
    
    ### based on moving threshold
    myVarThresh <- mqt(TS, Qdr, WinSize)
    
    for (i in 1:length(temp$Flow)) {
        temp$Thresh[i] <- myVarThresh[doy[i]]
        if (temp$Flow[i] < myVarThresh[doy[i]]) {
            temp$BelowThresh[i] <- TRUE}
        else {temp$BelowThresh[i]<-FALSE}
    } 
    
    return(temp)
}