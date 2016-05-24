##' Robust SNHT
##' 
##' This function performs a standard normal homogeneity test using a robust
##' estimator of the mean and standard deviation.  It also allows for a user-
##' defined definition of these statistics.
##' 
##' @usage robustSNHT(data, period, scaled=TRUE, rmSeasonalPeriod=Inf
##'       ,estimator=function(x, minObs=5){
##'       x = x[!is.na(x)]
##'       if(length(x)<minObs) #Too many NA values, don't return a result
##'         return(c(NA,NA))
##'       if(max(table(x))>length(x)/2) #Too many duplicate values, MAD will be 0
##'         return(c(NA,NA))
##'       fit = MASS::huber(x)
##'       return(c(fit[[1]], fit[[2]]))        
##' })
##' @param data The data to be analyzed for changepoints.
##' @param period The SNHT works by calculating the mean of the data on the
##' previous period observations and the following period observations. Thus,
##' this argument controls the window size for the test statistics.
##' @param estimator A custom function may be supplied to this function which
##' computes estimates for the mean and standard deviation. The function should
##' only take one argument (a numeric vector of data) and should return a
##' vector of length two: the estimated center and spread. The huber function
##' from MASS is implemented for the robust SNHT by default (along with some
##' data quality checks).
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
##' @importFrom zoo rollapply
##' 

robustSNHT <- function(data, period, scaled = TRUE, rmSeasonalPeriod = Inf
#      #Tukey estimator:
#      ,estimator=function(x){
#          fit = rlm( x ~ 1, psi=psi.bisquare, c=1.5, maxit=100)
#          return(c(fit$coeff, fit$s))        
#        }
      #Huber estimator:
      ,estimator = function(x, minObs = 5){
          #minObs arbitrarilty set to 5, want to ensure a decent number of values
          x = x[!is.na(x)]
          if(length(x) < minObs) #Too many NA values, don't return a result
            return(c(NA, NA))
          if(max(table(x)) > length(x) / 2) #Too many duplicate values, MAD will be 0
            return(c(NA, NA))
          fit = MASS::huber(x)
          return(c(fit[[1]], fit[[2]]))        
        }
      )  
{
  #Data quality checks
  if(period != round(period)){
    warning("period should be an integer!  Rounding to continue...")
    period = round(period)
  }
  if(!is.numeric(data))
    stop("data must be numeric!")
  if(2 * period >= length(data) - 1)
    stop("period is too large to compute statistics!")
  if(length(data) < 5)
    stop("snht requires at least 5 observations!")
  
  if(rmSeasonalPeriod < Inf)
    data = removeSeasonalPeriod(data, period = rmSeasonalPeriod)
  
#   #Difficult to implement appropriately, as NA's cause different number of values to occur
#   if(rmAC){
#     ACF = acf(d$data, na.action=na.pass, plot = F, lag.max = 2*period+1)$acf[-1]
#     corrAdj = 2*period + 
#       #N-1 elements with a lag of 1, N-2 elements with a lag of 2, etc.
#       4 * sum(ACF[1:(period-1)]*((period-1):1) ) -
#       #top row has lags N+1, ..., 2*N
#       #bottom row has lags 2, ..., N+1
#       #=> 1 lag of 2, 2 lags of 3, ..., N lags of N+1, N-1 lags of N+2, ..., 1 lag of 2*N
#       2*sum(ACF[2:(period+1)]*(1:period)) -
#       2*sum(ACF[(period+2):(2*period)]*((period-1):1))
#   } else {
#     corrAdj = 2*period
#   }
  
  #Compute rolling means for 1:period, 2:(period+1), etc.
  #by=2 calculates the stat for every day, i.e. every two obs.  Could speed up by increasing this value
  Means = zoo::rollapply(data, width=period, by=1, FUN=estimator)
  n = zoo::rollapply(data, width=period, by=1, FUN=function(x) sum(!is.na(x)) )
  #Right means start at observation period+2 (first obs to use is at period+1
  # and then means are to the right)
  rMeans = Means[(period + 2):nrow(Means), ]
  rN = n[(period + 2):nrow(Means)]
  #Left means are same length as right means but are at the beginning instead of end
  lMeans = Means[1:NROW(rMeans), ]
  lN = n[1:nrow(rMeans)]
  totMean = (lN*lMeans[, 1] + rN * rMeans[, 1]) / (lN + rN)
  if(scaled){
    scores = data.frame(
      #tukeyR also computes variances.  Variance of differences is sum of variances.
      #Original test stat has (N/2*(mu_L-mu)^2+N/2*(mu_R-mu)^2 )/sigma, but it assumes
      #N/2 observations on each side.  We must adjust for differing counts due to NA's:
      #Additional Note: Haimberger's test statistic divided by s but should have been s^2.
      #If scaled=TRUE, the correct statistic is used:
      score= (lN * (lMeans[, 1] - totMean)^2 + rN * (rMeans[, 1] - totMean)^2 ) /
        ( (lN * lMeans[, 2]^2 + rN * rMeans[, 2]^2) / (lN + rN) )
      ,leftMean = lMeans[, 1], rightMean = rMeans[, 1])
  } else {
    scores = data.frame(
    #tukeyR also computes variances.  Variance of differences is sum of variances.
    #Original test stat has (N/2*(mu_L-mu)^2+N/2*(mu_R-mu)^2 )/sigma, but it assumes
    #N/2 observations on each side.  We must adjust for differing counts due to NA's:
      score= (lN*(lMeans[,1]-totMean)^2 + rN*(rMeans[,1]-totMean)^2 ) /
        sqrt((lN*lMeans[,2]^2+rN*rMeans[,2]^2)/(lN+rN))
      ,leftMean=lMeans[,1], rightMean=rMeans[,1])
  }
  
  #Add zeros for rows skipped at beginning/end:
  toBind = data.frame(score=rep(NA,period), leftMean=NA, rightMean=NA)
  scores = rbind(toBind, scores, toBind)
  
  return(scores)
}
