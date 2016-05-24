##' Plot SNHT
##' 
##' Function to plot the result of the SNHT fit
##' 
##' @param data The vector of time series observations that were input to the
##' snht function.
##' @param stat A data.frame as returned by the snht function.
##' @param time If the observations in data are not equally spaced, then this
##' vector will specify the times of the observations.  This object should be
##' numeric or should be able to be coerced to numeric.
##' @param alpha The confidence level for the SNHT test.  Note, though, that
##' multiple tests are being performed and that is NOT accounted for in this
##' function.
##' 
##' @return No object is returned, but a plot is instead generated.
##' 
##' @export
##'
##' @import ggplot2
##' @importFrom gridExtra grid.arrange
##' @importFrom methods is
##' @importFrom stats qchisq
##' 

plotSNHT = function(data, stat, time = NULL, alpha = NULL){
    
    ## Data Quality Checks
    stopifnot(is.numeric(data))
    stopifnot(is(stat, "data.frame"))
    stopifnot("score" %in% colnames(stat))
    if(!is.null(alpha))
        stopifnot(is.numeric(alpha))
    if(is.null(time)){
        time = 1:length(data)
        stopifnot(nrow(stat) == length(data))
        if("time" %in% colnames(stat))
            stop("The snht statistic wasn't computed on evenly spaced data. ",
                 "Please supply the times to this plotting function.")
        stat$time = time
    } else {
        stopifnot("time" %in% colnames(stat))
        time = as.numeric(time)
        if(any(is.na(time)))
            stop("Supplied times were coerced to numeric, and the result has ",
                 "missing values.  Please check your time vector.")
    }
    
    pData = qplot(x = time, y = data)
    pStat = qplot(x = stat$time, y = stat$score, geom = "line") +
        labs(x = "time", y = "SNHT Statistic")
    if(!is.null(alpha))
        pStat = pStat + geom_hline(yintercept = qchisq(1-alpha, df = 1),
                                   linetype = 4, color = "blue")
    pStat = pStat + geom_vline(xintercept = stat$time[which.max(stat$score)],
                               color = "red", linetype = 4)
    pData = pData + geom_vline(xintercept = stat$time[which.max(stat$score)],
                               color = "red", linetype = 4)
    print(gridExtra::grid.arrange(pData, pStat))
}