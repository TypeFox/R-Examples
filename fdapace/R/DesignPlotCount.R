# This function is used to create count matrix based on
# observed pairs of time points for the raw covariance.

######
# Input:
######
#  t:  n * 1 array contains time points for n subjects
#  obsGrid: 1 * N vector contains sorted unique time points from t
#  noDiagonal: TRUE: set diagonal count as 0
#              FALSE: don't set diagonal count as 0
#  isColorPlot: TRUE: the resulting matrix has 1 indicates there exists points for out1(i) and out1(j)
#               FALSE: the resulting matrix contains counts of points for out1(i) and out1(j)
######
# Output: 
######
#  res: N * N matrix contains count for each distinct pairs of 
#             time points

DesignPlotCount = function(t, obsGrid, noDiagonal, isColorPlot){
  N = length(obsGrid) # number of distinct observed time pts
  res = matrix(0, nrow = N, ncol = N)

  for(cur in t){
    curidx = match(cur, obsGrid)
    if(isColorPlot == FALSE){
      res[curidx, curidx] = 1
    } else {
      res[curidx, curidx] = res[curidx, curidx] + 1
    }
  }

  if(noDiagonal == TRUE){
    diag(res) = 0
  }

  return(res)
}

# searchID = function(cur, obsGrid){
  # ni = length(cur)
  # id = rep(0, ni)
  # for(i in 1:ni){
    # id[i] = which(obsGrid == cur[i])
  # }
  # return(id)
# }
