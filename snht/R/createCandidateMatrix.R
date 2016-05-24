##' Create candidate matrix
##' 
##' This function creates candidate, a matrix where the (i,j)th entry
##' corresponds to the number of changepoints in difference series for location
##' j that occured at time i. For example, suppose location i_1 was paired with
##' i_2, i_3, and i_4.  If the statistic for i_1-i_2 and i_1-i_4 exceeded the
##' threshold at time j, then candidate_{i,j} = 2.
##' 
##' @param data The data.frame containing the observations, restructured as in
##' pairwiseSNHT.  So, the first column should be time, and the other columns
##' should be named with the locations and contain the observed values at each
##' location.
##' @param statistics The time x (number of pairs) matrix of SNHT statistics
##' computed for each difference series.
##' @param pairs The list object whose ith element specifies the neighboring
##' locations to the ith location.
##' @param crit The critical value such that if the snht statistic is larger
##' than crit, a changepoint is assumed to have occured.  Defaults to 100, as
##' recommended in Haimberger (see references).
##' 
##' @return A matrix of dimension time x (number of locations).  The (i,j)
##' element of this matrix indicates the number of changepoints found in
##' difference series containing the jth location at time i.
##' 

createCandidateMatrix = function(data, statistics, pairs, crit){
  locations = colnames(data)[-1]
  candidate = matrix(0, nrow=nrow(data), ncol=length(locations))
  colnames(candidate) = locations
  for(j in 1:ncol(statistics)){
    name = colnames(statistics)[j]
    name = strsplit(name, "-")[[1]]
    delta = as.numeric(statistics[, j] > crit)
    delta[is.na(delta)] = 0
    ## The ith list element of pairs specifies which stations are adjacent to 
    ## station i.  Thus, when updating the candidate matrix, we should update 
    ## column "a" if the pair "a"-"b" is significant and "b" is in pairs["a"]. 
    ## Also, we should update column "b" if "a" is pairs["b"].
    if(name[2] %in% pairs[name[1]][[1]])
      candidate[, name[1]] = candidate[, name[1]] + delta
    if(name[1] %in% pairs[name[2]][[1]])
      candidate[, name[2]] = candidate[, name[2]] + delta
  }    
  return(candidate)
}