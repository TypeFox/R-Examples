##' Standard Normal Homogeneity Test
##' 
##' This function performs a standard normal homogeneity test on the data 
##' supplied. This test searches the data for potential changepoints.
##' 
##' @usage snht(data, period, robust = F, time = NULL, scaled = TRUE, 
##'   rmSeasonalPeriod = Inf, ...)
##'   
##' @param data The data to be analyzed for changepoints.
##' @param period The SNHT works by calculating the mean of the data on the 
##'   previous period observations and the following period observations. Thus, 
##'   this argument controls the window size for the test statistics.
##' @param robust Flag indicating whether or not robust estimators should be 
##'   used. If T, then Huber's robust estimator for the mean and variance will
##'   be used (see ?MASS::huber).
##' @param time Numeric vector specifying times for the observations. If not 
##'   supplied, it is assumed that each observation occurs on one time period.
##'   If supplied, then the algorithm will create a new dataset with the same
##'   number of observations for each time unit by adding missing values.
##' @param scaled In the Haimberger paper, a typo is reported in the test 
##'   statistic.  The denominator ought to be s^2 instead of s, as the 
##'   distribution of the test statistic then becomes chi-squared assuming
##'   errors are normal.  If scaled = TRUE, the corrected version is used.  In
##'   this case, the test statistics can be compared to a chi-squared.  However,
##'   if scaled = FALSE, there is no clear distribution to compare with.
##' @param rmSeasonalPeriod This algorithm will overestimate the standard error 
##'   (and hence incorrectly estimate the test statistic) if there is strong 
##'   seasonality in the data.  By setting rmSeasonalPeriod to some value, a GAM
##'   model will be built to capture that seasonality.  Once it is estimated, it
##'   is removed and the SNHT statistic is computed on the resulting dataset. 
##'   Setting this argument to Inf prevents any modeling.
##' @param ... Other parameters, see ?robustSNHT, ?robustSNHTunequal.
##'   
##' @details The SNHT works by calculating the mean of the data on the previous 
##'   period and on the following period. The test statistic at each observation
##'   is then computed as described in Haimberger (2007). Essentially, though,
##'   it just compares the means of these two periods and normalizes by the
##'   standard deviation.
##'   
##'   Note: if there are not enough observations both before and after the 
##'   current observation, no test is performed.
##'   
##'   Large values of the test statistic suggests the presence of a changepoint.
##'   Haimberger (see references) suggests values larger than 100 should be 
##'   considered changepoints.  However, this does not apply if scaled = TRUE.
##'   
##'   Observations which are less than period away from the start or end of the 
##'   dataset do not have valid SNHT statistics.  Thus, the statistic for these 
##'   observations is returned as NA.
##'   
##' @return Returns a data.frame, with columns score, leftMean, and rightMean. 
##'   Statistic is the SNHT test statistic described above, and leftMean 
##'   (rightMean) are the means to the left (right) of the current observation.
##'   
##'   Additionally, if time is supplied, then time is returned on the output 
##'   data.frame. Note that new (missing) observations were introduced to the 
##'   dataset to ensure the same number of observations occur per day.
##'   
##' @author Josh Browning (jbrownin@@mines.edu), Carina Schneider
##'   (carina.schneider@@uzh.ch)
##'   
##' @references L. Haimberger. Homogenization of radiosonde temperature time 
##'   series using innovation statistics. Journal of Climate, 20(7): 1377-1403, 
##'   2007.
##'   
##' @seealso \code{\link{huber}}
##' @family snht functions
##' @examples
##' data = rnorm(1000)
##' brk = sample(1000, size=1)
##' data[1:brk] = data[1:brk]-2
##' out = snht( data, period=50, robust=FALSE )
##' summary(out)
##' 
##' data = rnorm(1000)
##' time = 1:1000 + rnorm(1000)
##' brk = sample(1000, size=1)
##' data[1:brk] = data[1:brk]-2
##' out = snht( data, period=50, time=time, robust=FALSE )
##' summary(out)
##' 
##' @keywords ~snht ~homogeneity
##'   
##' @export
##' @importFrom stats sd
##' 

snht = function(data, period, robust = F, time = NULL, scaled = TRUE
      ,rmSeasonalPeriod = Inf, ...){
  if(period != round(period)){
    warning("period should be an integer!  Rounding to continue...")
    period = round(period)
  }
  if(!is.numeric(data))
    stop("data must be numeric!")
  if(2 * period >= length(data) - 1 & is.null(time))
    stop("period is too large to compute statistics!")
  if(!is.null(time)){
    if(length(time) != length(data))
      stop("If time is not NULL, it must be the same length as data!")
    if(2 * period >= max(time) - min(time))
      stop("period is too large to compute statistics!")
  }
  if(rmSeasonalPeriod < Inf & rmSeasonalPeriod > length(data)/2)
    stop("Seasonal period must be <= half of the number of observations")
  if(length(data) < 5)
    stop("snht requires at least 5 observations!")
  
  if(!is.null(time)){
    if(robust)
      out = robustSNHTunequal( data = data, period = period, time = time, estimator = NULL
        ,scaled = scaled, rmSeasonalPeriod = rmSeasonalPeriod )
    else
      out = robustSNHTunequal(data=data, period=period, time=time, scaled=scaled
                              ,estimator=function(x){ #Use the mean and sd for non-robust fit
                                  x = x[!is.na(x)]
                                  return( c(mean(x), sd(x)) )
                              }, rmSeasonalPeriod = rmSeasonalPeriod )
    return(out)
  }
  
  if(!robust){
    #Use the "Robust" SNHT function, but supply a non-robust function: mean and sd.
    out = robustSNHT( data, period=period, scaled=scaled, estimator=function(x){
      x = x[!is.na(x)]
      return( c(mean(x), sd(x)) )
    }, rmSeasonalPeriod = rmSeasonalPeriod )
  } else {
    out = robustSNHT(data, period, scaled = scaled
        ,rmSeasonalPeriod = rmSeasonalPeriod, ...)
  }
  colnames(out) = c("score", "leftMean", "rightMean")
  return(out)
}