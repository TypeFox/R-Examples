#' Moving quantile threshold
#' 
#' This function calculates the daily moving window quantile threshold for use in 
#' identifying the partial duration series of streamflow droughts. 
#' @param TS output from \code{\link{create.ts}} containing a data.frame of flow
#'   time series
#' @param Qdr Numeric value of the drought threshold quantile.  Default is 0.2.
#' @param WinSize Numeric value specifying the size of the moving window in
#'   days.  Default is 30.
#' @details The threshold is defined by a moving quantile, where daily threshold 
#'   values are based on the 80th percentile of the flow duration curve 
#'   (i.e. 0.2 quantile) from a 30-day moving window (Beyene et al. 2014). 
#'   With this method, every day of the year has a different threshold based on 
#'   the streamflow measured on the day, the 15 days before 
#'   the day, and the 15 days after the day.The size of the moving window can be 
#'   modified with the WinSize argument, and the percentile can be modified with the
#'   Qdr argument.
#' @return Returns a numeric vector containing the streamflow drought threshold
#'   in m3/s for each day of the year.
#' @author Jennifer Dierauer
#' @references Beyene, B.S., Van Loon, A.F., Van Lanen, H.A.J., Torfs, P.J.J.F., 2014. 
#'   Investigation of variable threshold level approaches for hydrological drought 
#'   identification. Hydrol. Earth Syst. Sci. Discuss. 11, 12765-12797. 
#'   http://dx.doi.org/10.5194/hessd-11-12765-2014.
#' @seealso See \code{\link{create.ts}} to format the input flow series.
#' 
#'   The following functions use this function: \code{\link{dr.pds}}, 
#'   \code{\link{dr.events}}, \code{\link{dr.seas}}
#' @export
#' @examples
#' data(cania.sub.ts)
#' res <- mqt(cania.sub.ts)
#' 
#' # subset one year of the flow series
#' flow.sub <- cania.sub.ts[cania.sub.ts$year == 1990,]
#' 
#' # plot the 1990 observed flows in dark blue and the daily drought threshold in red
#' plot(flow.sub$doy, flow.sub$Flow, ylab="Q (m3/s)", xlab="Day of Year",
#'  pch=19, col="darkblue", type="b")
#' points(res, pch=19, cex=0.7, col="red")

mqt <- function(TS, Qdr=0.2, WinSize=30){
    
    doy <- TS$doy
    output <- array(NA, 365)
    doys <- c((365-(0.5*WinSize)):365, 1:365, 1:(0.5*WinSize))
        
    for (i in 1:365) {
        selected <- doys[i:(i+WinSize)]
        temp <- TS[doy %in% selected,]
        ## added if statement to handle seasonal 
        if (length(temp$Flow) > WinSize) {
            output[i] <- stats::quantile(temp$Flow, Qdr, na.rm=TRUE)
        } else {output[i] <- 0}
    }
    output[366] <- (output[365] + output[1])/2
    attr(output, "names") <- c(1:366)
    return(output)
}





