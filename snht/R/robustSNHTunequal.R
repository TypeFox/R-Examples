##' Robust SNHT with Unequal Times
##' 
##' This function performs a standard normal homogeneity test, but allows for
##' unequally spaced observations in time.  
##' 
##' @usage robustSNHTunequal(data, period, time, estimator = NULL, scaled=TRUE
##'     ,rmSeasonalPeriod = Inf)
##' 
##' @param data The data to be analyzed for changepoints.
##' @param period The SNHT works by calculating the mean of the data on the
##' previous period observations and the following period observations. Thus,
##' this argument controls the window size for the test statistics.
##' @param time Numeric vector specifying times for the observations. If not
##' supplied, it is assumed that each observation occurs on one time period. If
##' supplied, then the algorithm will create a new dataset with the same number
##' of observations for each time unit by adding missing values.
##' @param estimator See ?robustSNHT
##' @param scaled See ?snht.
##' @param rmSeasonalPeriod See ?snht.
##' 
##' @details The SNHT works by calculating the mean of the data on the previous
##' period and on the following period. The test statistic at each observation
##' is then computed as described in Haimberger (2007). Essentially, though, it
##' just compares the means of these two periods and normalizes by the standard
##' deviation.
##' 
##' Note: if there are not enough observations both before and after the
##' current observation, no test is performed.
##' 
##' Large values of the test statistic suggests the presence of a changepoint.
##' Haimberger (see references) suggests values larger than 100 should be
##' considered changepoints.  However, this does not apply if scaled = TRUE.
##' 
##' @return Returns a data.frame, with columns score, leftMean, and rightMean,
##' and time. Statistic is the SNHT test statistic described above, and leftMean
##' (rightMean) are the means to the left (right) of the current observation.
##' 
##' Note that new (missing) observations were introduced to the dataset to
##' ensure the same number of observations occur per day.
##' 
##' @author Josh Browning (jbrownin@@mines.edu)
##' 
##' @references L. Haimberger. Homogenization of radiosonde temperature time
##' series using innovation statistics. Journal of Climate, 20(7): 1377-1403,
##' 2007.
##' 
##' @seealso \code{\link{huber}}
##' @family snht functions
##' 
##' @keywords ~snht ~homogeneity ~robust
##' 
##' @export
##'
##' @importFrom plyr ddply
##' @importFrom methods is
##' 

robustSNHTunequal <- function(data, period, time, estimator=NULL, scaled=TRUE
                      ,rmSeasonalPeriod = Inf){
  #Data quality checks
  if(period != round(period)){
    warning("period should be an integer!  Rounding to continue...")
    period = round(period)
  }
  if(!is.numeric(data))
    stop("data must be numeric!")
  if(2 * period >= length(data) - 1 & is.null(time))
    stop("period is too large to compute statistics!")
  if(is.null(time)){
    if(2 * period >= max(time) - min(time))
      stop("period is too large to compute statistics!")
  }
  if(!is.numeric(time))
    stop("time must be numeric!")
  if(!is.null(estimator) & !is(estimator,"function"))
    stop("estimator must be a function or must be NULL!")
  if(length(data) != length(time))
    stop("data and time must be of the same length!")
  if(any(is.na(time)))
    stop("time cannot have missing values!")
  if(length(data) < 5)
    stop("snht requires at least 5 observations!")
  if(any(time != floor(time))){
    warning("Only integer values of time are used!  Rounding down.")
    time = floor(time)
  }

  d = data.frame(data=data, time=time, realObs=1)
  d = d[order(d$time),]
  #Bind on rows to d so we can see times that have no data as well
  dTemp = rbind(d, data.frame(data=NA, time=min(time):max(time), realObs=0) )
  obsPerDay = plyr::ddply(dTemp, "time", function(df){sum(df$realObs)})
  #Compute maximum observations for any time unit (called "day", but arbitrary)
  maxObs = max(obsPerDay[,2])
  obsPerDay$toBind = maxObs - obsPerDay[,2]
  #Determine which times are missing, and repeat the appropriate amount.  Then, bind on.
  toBind = rep( obsPerDay$time, times=obsPerDay$toBind )
  if(length(toBind)>0)
    d = rbind(d, data.frame(data=NA, time=toBind, realObs=0) )
  d = d[order(d$time),]
  
  if(is.null(estimator)){
    out = robustSNHT(data = d[, 1], period = period * maxObs, scaled = scaled
              ,rmSeasonalPeriod = rmSeasonalPeriod * maxObs)
  } else {
    out = robustSNHT(data = d[, 1], period = period * maxObs, scaled = scaled,
                     estimator, rmSeasonalPeriod = rmSeasonalPeriod * maxObs)
  }
  out$time = d$time
  return(out[d$realObs==1,])
}
