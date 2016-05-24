#' Calculate mean annual minimum n-day flows
#' 
#' This function calculates the calculates the mean annual minimum n-day flow by
#' calendar year or by hydrologic year.  This function can also be used to find
#' the annual minimum series by setting n=1.
#' @param TS output from \code{\link{create.ts}} containing a data.frame of flow
#'   time series
#' @param n Numeric value for the number of days in the n-day flow period. Default is 7.
#' @param by Character string indicating whether to use hydrologic years or 
#'   calendar years.  Default is "hyear".  Other option is "year".
#' @return Returns a numeric vector containing the calculated MAM n-day flow for each 
#'   year in the input time series.  The "times" attribute provides the corresponding
#'   year for each calculated value.
#' @author Jennifer Dierauer
#' @seealso \code{\link{screen.metric}}
#' @export
#' @examples
#' data(cania.sub.ts)
#' 
#' # find the annual minimum series and plot 
#' res <- MAMn(cania.sub.ts, n=1)
#' res2 <- screen.metric(res, "Q (m3/s)")
#' 
#' # do the same with MAM 7-day flow instead of annual minimum
#' res <- MAMn(cania.sub.ts, n=7)
#' res2 <- screen.metric(res, "Q (m3/s)")

MAMn <- function(TS, n=7, by="hyear") {
    
    ifelse(by=="hyear", Year <- as.factor(TS$hyear), Year <- as.factor(TS$year))
    
    YearList <- as.character(unique(Year))
    out <- (1:length(unique(Year)))
    attr(out, "times") <- YearList
    for (i in 1:length(YearList)) {
        temp <- TS[Year %in% YearList[i],]
        out[i] <- max(temp$Flow, na.rm=TRUE)
        
        for (j in 1:(length(temp$Flow)-n)) {
            checkmin <- mean(temp$Flow[j:j+(n-1)], na.rm=TRUE)
            if (j > 1) {
                if (checkmin < out[i]) {out[i] <- checkmin}
            }
        }
    }
    return(out)
}