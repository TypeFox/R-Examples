###############################################################################
### Description: Identifies all paths length less than or equal to max.length
###              between all pairs of competitors.  Used in conductance.
### Input:
###   conf - N-by-N conflict matrix whose (i,j)th element is the number of 
###          times i defeated j
###   maxLength - positive numeric integer indicating the maximum length
###                of paths to identify
### Output: list whose elements are all paths of a given length
###############################################################################

allPaths = function(conf, maxLength){
  paths = lapply(2:maxLength, FUN = function(l, conf){
    do.call(rbind, lapply(1:nrow(conf), FUN = IDpaths, conf = conf, l = l))
  }, conf = conf)  
  return(list(which(conf > 0, arr.ind = TRUE), paths))
}