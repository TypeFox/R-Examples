##' Get Unique Pairs
##' 
##' For the pairwise SNHT, many difference time-series must be computed.
##' However, if location i is used for location j, then it's very likely that
##' location j will be used for location i.  Thus, the pairs object will likely
##' have many duplicate pairs.  To save computation time, this function finds
##' which pairs are unique.
##' 
##' @param pairs The pairs list object, as returned by ?getPairs.
##' 
##' @return data.frame with columns loc1 and loc2.  This data.frame will have
##' no duplicates, and describes all the pairs that need to be computed.
##' 

getUniquePairs = function(pairs){
  #Compute all pairs that will need to be analyzed:
  uniquePairs = data.frame(loc1 = rep(names(pairs), lapply(pairs,length))
                          ,loc2 = do.call("c", pairs), stringsAsFactors=F )
  #Remove non-unique pairs.  Start at nrow(pairsDf) and work down.  Otherwise, the
  #iterator would break when you delete a row, as the nrow(pairsDf) would reduce but
  #the iterator would still run to the original nrow(pairsDf)
  for(i in nrow(uniquePairs):1)
    if(any(uniquePairs[,1]==uniquePairs[i,2] & uniquePairs[,2]==uniquePairs[i,1]))
      uniquePairs = uniquePairs[-i,]
  return(uniquePairs)
}