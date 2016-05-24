#' Seasonal baseflow percentage
#' 
#' This function estimates the percentage of baseflow in a given period relative to the total
#' annual baseflow.
#' @param TS output from \code{\link{create.ts}} containing a data.frame of flow
#'  time series
#' @param seas Integers representing months of the year. Default is c(6:8), i.e. June-August.
#' @details This function calls \code{\link{bf_eckhardt}} to complete the 
#'   baseflow separation.
#' @return Returns a vector containing the calculated percentage for each year
#'   in the input time series.  The "times" attribute provides the corresponding
#'   year for each calculated value.
#' @author Jennifer Dierauer
#' @seealso See \code{\link{bf.stats}} to calculate additional baseflow metrics.
#' @export
#' @examples
#' data(cania.sub.ts)
#' res <- bf.seas(cania.sub.ts)
#' res2 <- screen.metric(res, "Percent Annual Baseflow in Jun-Aug")

bf.seas <- function(TS, seas=c(6:8)) {
    
    ### Set parameter values for Eckhardt RDF
    BFindex <- 0.8
    alpha <- 0.970 ##based on values suggested by Eckhardt 2012 for perennial stream
    
    ## calculate daily BF
    TS <- subset(TS, !is.na(TS$Flow))
    TS$bf <- bf_eckhardt(TS$Flow, alpha, BFindex) 
    
    ## loop through years and calculate seasonal baseflow percentage
    years <- unique(TS$year)
    BFpct <- numeric()
    for (i in 1:length(years)) {
        
        TS.sub <- subset(TS, TS$year == years[i])
        BFvol <- sum(TS.sub$bf)
        TS.sub <- subset(TS, TS$month %in% seas)
        
        if (length(TS.sub) == 1) {
            BFpct.sub <- 0
        } else {
            BFpct.sub <- (sum(TS.sub$bf)/BFvol)
        }
        
        BFpct <- c(BFpct, BFpct.sub)
    }
    
    attr(BFpct, "times") <- as.numeric(years)
    attr(BFpct, "dimnames") <- NULL
    
    return(BFpct)
    
}