#' Sum missing data points from a daily time series
#' 
#' Counts the number of missing data points by calendar year, hydrologic year, or month
#' @param input output from \code{\link{NA.runs}}
#' @param by character string identifying the time period to summarize by.  
#'   Defaults is hydrologic year ("hyear"), other choices are "year" and "month".
#'   The "month" option will return the number of missing data points for each
#'   month in the time series.
#' @param hyrstart optional argument, define start month of hydrologic year
#' @return Returns a numeric vector of the number of missing observations per
#'   summary period.  The "times" attribute of the returned vector provides the
#'   corresponding year, hyear, or month.
#' @author Jennifer Dierauer
#' @seealso \code{\link{NA.runs}}
#' @export
#' @examples
#' data(caniapiscau)
#' cania.sub <- caniapiscau[300:1200,]
#' cania.ts <- create.ts(cania.sub)
#' res <- NA.runs(cania.ts)
#' res2 <- NA.sum(res)


NA.sum <- function(input, by="hyear", hyrstart=1) {
    
    NAruns <- as.Date(NA)
    for (i in 1:length(input$Start)) {
        mseq <- seq(from=input$Start[i], to=input$End[i], by=1)
        NAruns <- c(NAruns, mseq)
    }
    NAruns <- NAruns[-1]
    NAruns <- data.frame(Date=NAruns)
    NAruns <- YMD.internal(NAruns)
    
    if (by == "year") {
        mcount <- tapply(NAruns$Date, NAruns$year, length)
        labs <- unique(NAruns$year)
        attr(mcount, "times") <- labs
    }
    
    if (by == "hyear") {
        NAruns <- hyear.internal(NAruns, hyrstart)
        mcount <- tapply(NAruns$Date, NAruns$hyear, length)
        labs <- unique(NAruns$hyear)
        attr(mcount, "times") <- labs
    }
    
    if (by == "month") {
        NAruns$yrm <- paste(NAruns$year, NAruns$month, sep="-")
        mcount <- tapply(NAruns$Date, NAruns$yrm, length)
        labs <- unique(NAruns$yrm)
        attr(mcount, "times") <- labs
    }
    
    return(mcount)
}


