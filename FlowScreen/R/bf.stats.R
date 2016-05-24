#' Baseflow statistics
#' 
#' This function estimates the baseflow and calculates the mean, max, and min
#' baseflow and baseflow index for a user defined time period.
#' @param TS output from \code{\link{create.ts}} containing a data.frame of flow
#'   time series
#' @param by summary period.  Options are "year", "hyear", "month", or "doy".  
#'   Default is "hyear".
#' @details This function calls \code{\link{bf_eckhardt}} to complete the 
#'   baseflow separation.
#' @return Returns a data.frame with the following columns:
#'   \itemize{
#'     \item By - Unique values representing the summary periods, e.g. a list of
#'       unique years, months, or days of year 
#'     \item MeanQ - Mean daily streamflow for the summary period, in m3/s 
#'     \item MeanBF - Mean daily baseflow for the summary period, in m3/s 
#'     \item MaxBF - Maximum daily baseflow for the summary period, in m3/s 
#'     \item MinBF - Minimum daily baseflow for the summary period, in m3/s 
#'     \item BFVol - Baseflow volume for the summary period, in km3 
#'     \item MeanBFI - Mean daily baseflow index for the summary period, dimensionless 
#'     \item MaxBFI - Maximum daily baseflow index for the summary period, dimensionless 
#'     \item MinBFI - Minimum daily baseflow index for the summary period, dimensionless 
#'   }
#' @author Jennifer Dierauer
#' @export
#' @examples
#' data(cania.sub.ts)
#' 
#' res <- bf.stats(cania.sub.ts)
#' res2 <- screen.metric(res[,2], "m3/s")

bf.stats <- function(TS, by="hyear") {
    
    ### Set parameter values for Eckhardt RDF
    BFindex <- 0.8
    alpha <- 0.970 ##based on values suggested by Eckhardt 2012 for perennial stream
    
    ## calculate daily BF and BFI
    ##BF filter can't handle NAs
    ##temporary fix
    
    temp <- subset(TS, !is.na(TS$Flow))
    temp$bf <- bf_eckhardt(temp$Flow, alpha, BFindex)  
    temp$bfi <- temp$bf/temp$Flow
    
    #create factors
    if (by=="hyear") {SumBy <- as.factor(temp$hyear)}
    if (by=="year") {SumBy <- as.factor(temp$year)}
    if (by=="month") {SumBy <- as.factor(temp$month)}
    if (by=="doy") {SumBy <- as.factor(temp$doy)}
    
    MyTitles <- c("MeanQ", "MeanBF", "MaxBF", "MinBF",
                  "BFVol", "MeanBFI", "MaxBFI", "MinBFI")
    
    ## calculate Mean Q, Max BF, BF Volume, and BFI
    output<-data.frame(By=as.numeric(as.character(unique(SumBy))))
    output[,2] <- tapply(temp$Flow, SumBy, mean) ## Mean Q
    output[,3] <- tapply(temp$bf, SumBy, mean) ## Mean BF
    output[,4] <- tapply(temp$bf, SumBy, max) ## Max Baseflow
    output[,5] <- tapply(temp$bf, SumBy, min) ## Min BF
    output[,6] <- tapply(temp$bf, SumBy, sum)
    output[,6] <- output[,4] * 0.031536000 ## Baseflow Volume (km3)
    output[,7] <- tapply(temp$bfi, SumBy, mean, na.rm=T) ## Mean BFI
    output[,8] <- tapply(temp$bfi, SumBy, max, na.rm=T)  ## Max BFI
    output[,9] <- tapply(temp$bfi, SumBy, min, na.rm=T)  ## Min BFI
    
    colnames(output)<- c(by, MyTitles)
    
    years <- as.character(output[,1])
    for (i in 2:9) {
        attr(output[,i], "times") <- years
    }
    
    return(output)
}