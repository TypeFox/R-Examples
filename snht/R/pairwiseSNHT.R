##' Pairwise Standard Normal Homogeneity Test
##' 
##' This function performs a pairwise standard normal homogeneity test on the
##' data supplied, as described in Menne & Williams (2009).
##' 
##' @usage pairwiseSNHT(data, dist, k, period, crit=100, returnStat=FALSE, ...)
##' @param data The data to be analyzed for changepoints.  It must be a
##' data.frame and contain either two or three columns.  The mandatory columns
##' are data and location, named as such.  The option column is time, and this
##' argument will be passed to snht.
##' @param dist A distance matrix which provides the distance between location
##' i and location j.  Rows and columns must be named with the locations in
##' data.  Note that non-symmetric distances may be used.  In that case,
##' neighbors for station i will be determined by the smallest values in the
##' row of dist corresponding to i.
##' @param k How many of the nearest neighbors should be used to construct
##' pairwise difference time series?  Note that more than k neighbors may be
##' used if there are ties in the distances between locations.
##' @param period The SNHT works by calculating the mean of the data on the
##' previous period observations and the following period observations.  Thus,
##' this argument controls the window size for the test statistics.
##' @param crit The critical value such that if the snht statistic is larger
##' than crit, a changepoint is assumed to have occured.  Defaults to 100, as
##' recommended in Haimberger (see references).
##' @param returnStat See return value.  If TRUE, the snht statistics for each
##' time point and for each difference pair are returned.
##' @param ... Additional arguments to pass to the snht function (such as
##' robust, time, or estimator).
##' 
##' @details The pairwise snht works with a set of time series.  For each time
##' series, it's closest k neighbors are determined, and a time series of the
##' difference between each of those time series is created.  The snht is then
##' applied to each of these difference time series.  Changepoints in one time
##' series can be detected by searching for large values of the test statistic
##' across all difference time series for a particular location.
##' 
##' The usefulness of the pairwise snht is that it removes any patterns in the
##' data that could affect the basic snht.  For example, seasonal and linear
##' trends that exist globally will be removed from the difference series, and
##' thus changepoints are more easily detected.
##' 
##' @return If returnStat is TRUE, the snht statistics for each time point and
##' for each difference pair are returned.
##' 
##' Otherwise, a named list is returned.  The first element, data, contains the
##' homogenized data in the same format as the supplied data.  The second
##' element, breaks, contains a data.frame where the first column is the
##' location where a break occured, the second column is the time of the break,
##' and the third column is the amount that data after the break was shifted
##' by.
##' 
##' @references L. Haimberger. Homogenization of radiosonde temperature time
##' series using innovation statistics. Journal of Climate, 20(7): 1377-1403,
##' 2007.
##' 
##' Menne, M. J., & Williams Jr, C. N. (2009). Homogenization of temperature
##' series via pairwise comparisons. Journal of Climate, 22(7), 1700-1717.
##' 
##' @author Josh Browning (jbrownin@@mines.edu)
##' keyword ~snht ~homogeneity ~pairwise
##' 
##' @export
##'
##' @importFrom plyr ddply
##' @importFrom reshape2 dcast
##' @importFrom reshape2 melt
##' @importFrom methods is
##' 

pairwiseSNHT <- function(data, dist, k, period, crit=100, returnStat=FALSE,
    ...){
  #data quality checks
  stopifnot(is(data,"data.frame"))
  if(is(data, "data.table")){
    stop("data must be a data.frame, not a data.table")
  }
  if(ncol(data)==2){
    stopifnot(colnames(data) %in% c("data","location"))
    # Reorder columns
    data = data[, c("data", "location")]
  }
  if(ncol(data)==3){
    stopifnot(colnames(data) %in% c("data","location","time"))
    # Reorder columns
    data = data[, c("data", "location", "time")]
  }
  stopifnot(ncol(data) %in% c(2,3))
  locs = as.character(unique(data$location))
  stopifnot(rownames(dist) == colnames(dist))
  stopifnot(all(rownames(dist) %in% locs))
  stopifnot(all(locs %in% rownames(dist)))
  stopifnot(k >= 1) #Must have at least one neighbor
  stopifnot(k <= length(locs)-1) #Can have at most length(locs)-1 neighbor, since self can't be used  
  stopifnot(diag(dist) == 0)
  if(any(dist[row(dist) != col(dist)]<=0))
    stop("Off diagonal elements of dist must be >0")
  
  pairs = getPairs(dist, k=k)
  uniquePairs = getUniquePairs(pairs)
  
  #Add times if they don't already exist (just 1:nrow()).
  if(!"time" %in% colnames(data)){
    if(length( unique( table( data$location ) ) ) != 1){
      stop("All locations must have the same number of obs if time is not provided!
           May need to remove unused levels in data.")
    }
    data$order = 1:nrow(data) #ensure original ordering is preserved
    data = plyr::ddply(data, "location", function(df){
      df = df[order(df$order),]
      df$time = 1:nrow(df)
      return(df)
    } )
    data$order = NULL
  }
  
  #Restructure data
  data = reshape2::dcast(data, formula = time ~ location, value.var = "data")
  diffs = data.frame(time = data$time)
  for(i in 1:nrow(uniquePairs)){
    diffs = cbind(diffs, data[,uniquePairs[i,1]] - data[,uniquePairs[i,2]])
    colnames(diffs)[ncol(diffs)] = paste0(uniquePairs[i,1],"-",uniquePairs[i,2])
  }
  
  #Compute snht statistics
  statistics = apply(diffs[,-1], 2, snht, period=period, time=diffs[,1], ...)
  avgDiff = do.call("cbind", lapply(statistics, function(x) x$rightMean-x$leftMean ) )
  statistics = do.call("cbind", lapply(statistics, function(x) x$score))
  if(returnStat)
    return(statistics)
  
  candidate = createCandidateMatrix(data, statistics = statistics,
                                    pairs = pairs, crit = crit)
  out = unconfoundCandidateMatrix(candidate = candidate, pairs = pairs,
    statistics = statistics, data = data, period = period, avgDiff = avgDiff)
  
  out$data = reshape2::melt(data = out$data, id.vars = "time")
  rownames(out$data) = NULL
  colnames(out$data)[colnames(out$data) == "value"] = "data"
  colnames(out$data)[colnames(out$data) == "variable"] = "location"
  out$data = out$data[,c("data", "location", "time")]
  return(out)
}