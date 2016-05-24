#' Day of year for annual maximum series
#' 
#' This function returns the day of the hydrologic year for each annual maximum flow.
#' @param TS output from \code{\link{create.ts}} containing a data.frame of flow
#'   time series
#' @return Returns a numeric vector containing the day of the (hydrologic) year for 
#'   each annual maximum flow. The "times" attribute contains the 
#'   hydrologic year for each element in the vector.
#' @author Jennifer Dierauer
#' @seealso See \code{\link{create.ts}} to format the input flow series.
#' 
#'   See \code{\link{pk.max}} for the annual maximum flow series.
#' @export
#' @examples
#' data(cania.sub.ts)
#' res <- pk.max.doy(cania.sub.ts)
#' res2 <- screen.metric(res, "Day of Year")


pk.max.doy <- function(TS) {
    
    year <- as.factor(TS$hyear)
    year_list <- unique(year)
    max.doy <- rep(NA, length(year_list))
    
    for (i in 1:length(year_list)) {
        
        temp <- subset(TS, TS$hyear %in% year_list[i])
        temp <- subset(temp, temp$Flow == max(temp$Flow))
        
        # if more than one day with same max flow value, take first
        if (length(temp$Flow > 1)) {temp <- temp[1,]}
        max.doy[i] <- as.numeric(temp$hdoy)
        
    }
    
    attr(max.doy, "times") <- as.numeric(as.character(year_list))
    return(max.doy)
}