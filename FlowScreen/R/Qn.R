#' Calculate flow quantiles
#' 
#' This function calculates flow quantiles by hydrologic year, calendar
#' year, month, or doy.
#' @param TS output from \code{\link{create.ts}} containing a data.frame of flow
#'   time series
#' @param n Numeric value of the quantile.  Default is 0.1.
#' @param by Character string indicating time unit to summarize by.  Default is
#'   "hyear" for hydrologic year, see \code{\link{create.ts}}.  Other options 
#'   are "year" for calendar year, "month", or "doy" for day of year. 
#' @return Returns a numeric vector of the calculated flow quantile for the time
#'   periods indicated with the "by" argument.  The "times" attribute contains the
#'   hydrologic year, calendar year, month, or day of year for each data point.
#' @author Jennifer Dierauer
#' @export
#' @examples
#' data(cania.sub.ts)
#' 
#' # 50% quantile, i.e. mean, by calendar year
#' res <- Qn(cania.sub.ts, n=0.5, by="year")
#' res2 <- screen.metric(res, "Q (m3/s)")
#' 
#' # Default 10% quantile, by hydrologic year
#' res <- Qn(cania.sub.ts)
#' res2 <- screen.metric(res, "Q (m3/s)")


Qn <- function(TS, n=0.1, by="hyear") {
  
    if (by=="hyear"){
        SumBy <- as.factor(TS$hyear)
        mList <- as.character(unique(SumBy))
    }
    
    if (by=="year"){
        SumBy <- as.factor(TS$year)
        mList <- as.character(unique(SumBy))
    }
    
    if (by=="month"){
        SumBy <- as.factor(TS$month)
        mList <- as.character(unique(SumBy))
    }
    
    if (by=="doy"){
        SumBy <- as.factor(TS$doy)
        mList <- as.character(unique(SumBy))
    }
    
    out <- tapply(TS$Flow, SumBy, stats::quantile, n)
    
    attr(out, "times") <- mList
    attr(out, "dimnames") <- NULL
    return(out)
}