#' Calculate the inter-event duration
#' 
#' This function calculates duration (in days) between flow peaks.
#' @param Peaks Output from \code{\link{pks}}.
#' @author Jennifer Dierauer
#' @return Returns a numeric vector containing the duration (in days) between peaks
#'   over threshold from \code{\link{pks}}. The "times" attribute contains the calendar 
#'   year date of the earlier peak. The "names" attribute contains the hydrologic year and
#'   the day (1-366) of the hydrologic year.
#' @export
#' @examples
#' data(cania.sub.ts)
#' res1 <- pks(cania.sub.ts)
#' res2 <- pks.dur(res1)
#' res3 <- screen.metric(res2, "Inter-Event Duration (days)")

pks.dur <- function(Peaks) {
    
    MyDates <- as.Date(attr(Peaks,"times")) # get dates from atomic vector
    hyr.hdoy <- attr(Peaks, "names")
    
    Durations <- c(2:length(MyDates)) # create vector to fill
    
    ## calculate inter-event duration
    for (i in 1:length(Durations)) {
        Durations[i] <- (MyDates[i+1] - MyDates[i])
    }
    
    ## attach time attribute
    attr(Durations, "times") <- MyDates[-1]
    attr(Durations, "names") <- hyr.hdoy[-1]
    
    return(Durations)
}