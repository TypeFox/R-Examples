#' Create a Time Series of daily streamflow observations
#' 
#' This function creates a daily time series formatted
#' for use with the functions in this package.
#' @param Flows Data.frame containing daily streamflow time series loaded
#'   with the \code{\link{read.flows}} function.
#' @param hyrstart define start month of hydrologic year. Defaults to 10 (October).
#' @return Returns a data.frame with year, month, doy, and hyear columns 
#'   appended to the original input data.frame.
#' @author Jennifer Dierauer
#' @export 
#' @examples
#' data(caniapiscau)
#' # subset flow series for shorter example run time
#' caniapiscau.sub <- caniapiscau[300:1800,]
#' caniapiscau.sub.ts <- create.ts(caniapiscau.sub)


create.ts <- function(Flows, hyrstart=10) {
    
    numobs <- length(Flows$Flow)
    mseq <- seq(from=Flows$Date[1], to=Flows$Date[numobs], by=1)
    mseq <- data.frame(Date=mseq, Flow=NA)
    
    if (length(mseq) > numobs) {
        out <- merge(mseq, Flows, by = "Date", all.x=T, all.y=T)
        out <- data.frame(ID=out$ID, Date=out$Date, Flow=out$Flow.y, 
                          Code=out$SYM, Agency=out$Agency)
        
    } else {out <- data.frame(ID=Flows$ID, Date=Flows$Date, Flow=Flows$Flow, 
                              Code=Flows$SYM, Agency=Flows$Agency)}
    
    out <- YMD.internal(out)
    out <- hyear.internal(out, hyrstart)
    
    out <- subset(out, !is.na(out$Flow))
    
    return(out)   
}