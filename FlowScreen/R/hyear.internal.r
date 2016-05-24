#' Add hydrologic Year, month, and doy columns to a daily time series
#' 
#' @param TS Output from \code{\link{create.ts}} function.
#' @param hyrstart define start month of hydrologic year. Defaults to 10 (October).
#' @return Returns a data.frame with hyear, hmonth, and hdoy columns 
#'   appended to the original input data.frame.
#' @author Jennifer Dierauer


hyear.internal <- function(TS, hyrstart=1) {
    
    TS$hyear <- TS$year
    TS$hmonth <- TS$month
    TS$hdoy <- TS$doy
    
    if (hyrstart==1) {TS$hyear <- TS$year}
    
    if (hyrstart > 6.5) {
        
        MonthsUp <- c(hyrstart:12)
        month.hyr <- data.frame(hyr=c(1:12), original=c(hyrstart:12, 1:(hyrstart-1)))
        
        for (i in 1:length(TS$Date)) {
            if (TS$month[i] %in% MonthsUp) {
                TS$hyear[i] <- TS$year[i] + 1
                hyr.start.date <- as.Date(paste(TS$year[i], hyrstart, "01", sep="-"),
                                          format="%Y-%m-%d")
            } else {
                hyr.start.date <- as.Date(paste((TS$year[i]-1), hyrstart, "01", sep="-"),
                                          format="%Y-%m-%d")
            }
            
            TS$hmonth[i] <- subset(month.hyr, month.hyr$original == TS$month[i])$hyr
            TS$hdoy[i] <- as.numeric(TS$Date[i] - hyr.start.date) + 1
        }
    }
    
    if (hyrstart < 6.5 & hyrstart != 1) {
        
        MonthsDown <- c(1:(hyrstart-1))
        month.hyr <- data.frame(hyr=c(1:12), original=c(hyrstart:12, 1:(hyrstart-1)))
        
        for (i in 1:length(TS$Date)) {
            if (TS$month[i] %in% MonthsDown) {
                TS$hyear[i] <- TS$year[i] - 1
                hyr.start.date <- as.Date(paste((TS$year[i] - 1), hyrstart, "01", sep="-"),
                                          format="%Y-%m-%d")
            } else {
                hyr.start.date <- as.Date(paste(TS$year[i], hyrstart, "01", sep="-"),
                                          format="%Y-%m-%d")
            }
            
            TS$hmonth[i] <- subset(month.hyr, month.hyr$original == TS$month[i])$hyr
            TS$hdoy[i] <- as.numeric(TS$Date[i] - hyr.start.date) + 1
        }
    }
    
    return(TS)
}
