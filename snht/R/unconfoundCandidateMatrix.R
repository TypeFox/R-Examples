##' Unconfound candidate matrix
##' 
##' This function "unconfounds" the candidate matrix.  At each time point and 
##' for each location, we have the number of difference series which resulted in
##' a changepoint.  The location with the largest count is assumed to be the
##' location where the changepoint occurs.  Assignment of changepoints should
##' then proceed iteratively, where each new changepoint is assigned based on
##' the current highest count.
##' 
##' @param candidate The candidate matrix, as computed by 
##'   ?createCandidateMatrix.
##' @param pairs The list object whose ith element specifies the neighboring 
##'   locations to the ith location.
##' @param statistics The time x (number of pairs) matrix of SNHT statistics 
##'   computed for each difference series.
##' @param data The data.frame containing the observations, restructured as in 
##'   pairwiseSNHT.  So, the first column should be time, and the other columns 
##'   should be named with the locations and contain the observed values at each
##'   location.
##' @param period The SNHT works by calculating the mean of the data on the 
##'   previous period observations and the following period observations.  Thus,
##'   this argument controls the window size for the test statistics.
##' @param avgDiff A matrix containing the average differences between time 
##'   series pairs.  Generally this is created within pairwiseSNHT().
##'   
##' @return A list of two elements.  The first element contains the data after 
##'   the breaks have been removed.  The second element is a data.frame with 
##'   information regarding the detected changepoints (or NULL if none are
##'   found).
##'   

unconfoundCandidateMatrix = function(candidate, pairs, statistics, data,
                                     period, avgDiff){
  #Algorithm:
  #maxCols = locations which all attain the max count
  #If (|maxCols|>1)
  #  Examine all difference series for all elements of maxCols but restricted to rows
  #  where the max count occurs.  Find the statistic with the largest value.  Call the
  #  two locations forming this series col_A, col_B.  Set the break time to the
  #  time of this statistic, call it brkT.
  #    brkCol = col_B or col_A, respectively
  #  else
  #    At time brkT, examine the other difference statistics for col_A and col_B.  Assign
  #    a break to whichever pair has the second largest statistic at time brkT
  #else
  #  brkCol = column with largest count
  #  Amongst all points with max count, find the time where the largest statistic
  #  occurs.  Call this time brkT
  #
  #Update candidate: Set t_window=brkT + -period:period.
  #candidate[t_window, brkCol] = 0
  #brkCol_pairs = set of all time series that use maxCol in their difference series.
  #candidate[t_window, brkCol_pairs] = pmax( candidate[t_window, brkCol_pairs]-1, 0)
  #Note: statistics must also be updated, to ensure that times from other stations are
  #not choosen when the difference is due to the homogenized station.
  #The assumption is that the changepoint in the maxCol series was the issue in the 
  #other difference series.  This seems pretty reasonable.
  
  breaks = NULL
  while(any(candidate>0)){
    colMax = apply(candidate, 2, max)
    maxVal = max(colMax)
    maxCols = colMax==maxVal
    maxCols = names(maxCols)[maxCols]
    
    #Determine which difference series we'll need to examine
    pairsMax = pairs[maxCols]
    pairsMax = data.frame(loc1 = rep(names(pairsMax), lapply(pairsMax,length))
                         ,loc2 = do.call("c", pairsMax) )
    #Could occur in either order, create both possibilities
    pairsMax$diff = paste0(pairsMax[,1], "-", pairsMax[,2])
    pairsMax$diff = ifelse(!pairsMax$diff %in% colnames(statistics)
                          ,paste0(pairsMax[,2], "-", pairsMax[,1])
                          ,pairsMax$diff )
    
    #Determine brkT, the time at which the current break is detected, and brkCol.
    #Note: if columns A and B have the same value in candidate and the maximum is
    #between A and B, then the first column will be used.  Tie-breaking in this case
    #is non-trivial (what if each only have candidate=1?  how can we pull other values?)
    brkTimes = lapply(maxCols, function(location){
      if(!location %in% maxCols)
        return(NULL)
      rows = which(candidate[, colnames(candidate) == location] == maxVal)
      cols = pairsMax$diff[pairsMax$loc1 == location]
      cols = colnames(statistics) %in% cols
      currCol = which.max(apply(statistics[rows,cols,drop=F], 2, max))
      currColTime = rows[which.max(statistics[rows,currCol])]
      return(data.frame(time = currColTime, stat = max(statistics[currColTime,])))
    })
    brkTimes = do.call("rbind", brkTimes)
    brkT = brkTimes$time[which.max(brkTimes$stat)]
    brkCol = maxCols[which.max(brkTimes$stat)]
    
    #Update candidate matrix
    tWindow = brkT + -period:period #No changepoints within period observations
    candidate[tWindow, brkCol] = 0
    adjCandCols = lapply(pairs, function(x){brkCol %in% x})
    adjCandCols = names(pairs)[do.call("c",adjCandCols)]
    candidate[tWindow, adjCandCols] = pmax(candidate[tWindow, adjCandCols]-1, 0)
    
    #Update statistics matrix
    adjStatCols = c(paste0(brkCol,"-",adjCandCols), paste0(adjCandCols,"-",brkCol))
    adjStatCols = adjStatCols[adjStatCols %in% colnames(statistics)]
    #Set to zero as we're assuming large values are caused by brkCol
    statistics[tWindow, adjStatCols] = 0
    
    #Update data
    adjMeanCols = c(paste0(brkCol,"-",pairs[[brkCol]]), paste0(pairs[[brkCol]],"-",brkCol))
    adjMeanCols = adjMeanCols[adjMeanCols %in% colnames(statistics)]
    currentDiff = avgDiff[brkT,adjMeanCols]
    #If brkCol isn't first in a pair, the value should be negated (as it's a
    #difference between brkCol and the pair station that we're interested in):
    reverseCols = grep(paste0(brkCol, "$"), names(currentDiff))
    currentDiff[reverseCols] = -currentDiff[reverseCols]
    shift = mean(currentDiff, na.rm=T)
    data[brkT:nrow(data),brkCol] = data[brkT:nrow(data),brkCol] - shift
    
    #Append detected break to breaks
    breaks = rbind(breaks, data.frame(brkT, brkCol, shift))
  }
  if(!is.null(breaks)){
    breaks = data.frame(breaks)
    colnames(breaks) = c("time", "location", "shift")
    breaks$time = as.numeric(breaks$time)
    breaks$shift = as.numeric(breaks$shift)
  }
  
  return(list(data = data, breaks = breaks))
}