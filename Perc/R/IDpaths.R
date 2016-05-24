### Need the function IDpaths() to find all dominance paths.

IDpaths = function(conf, i, len){
  ### IDpaths: function that identifies all unique dominance paths of order (len - 1) beginning at subject i
  ### conf: N-by-N conflict matrix whose (i,j)th element is the number of times i defeated j
  ### i: the subject at the beginning of each dominance path
  ### len: the length of the dominance paths to be identified (len = order + 1)
  
  if(sum(conf[i,] > 0) == 0) return(matrix(0, 0, len+1))
  # i = 1                            ###? does i always == 1?
  #len = 5
  levels = list()
  levels[[1]] = which(conf[i,] > 0)  ###for the ith individual, index the individual column location where this individual win over this individual directly.
  
  # find column index of individuals whom the ith individual win over through pathway of len.
  
  for(j in 2:len){
    levels[[j]] = lapply(unlist(levels[[j-1]]), function(k) which(conf[k,] > 0))
  }
  
  # create a matrix ret to represent the pathway information
  ret = matrix(0, nrow = length(unlist(levels[[len]])), ncol = len+1)
  
  # store information
  # the first column will always be the ith individual
  ret[,1] = i
  
  # the len+1th column will be individuals in levels[[len]]
  ret[,len+1] = unlist(levels[[len]])
  
  # if len == 2, the 2nd column will be items in levels[[1]] which
  #  find pathways in levels[[2]].
  if(len == 2){
    ret[,2] = rep(unlist(levels[[1]]), sapply(levels[[2]], length))
  }
  for(j in len:2){
    #j = 4
    currLengths = sapply(levels[[j]], length)
    if(j < len){
      effLengths = numeric(length(currLengths))
      ctr = 1
      for(d in 1:length(effLengths)){
        if(currLengths[d] != 0){
          effLengths[d] = sum(prevLengths[ctr:(ctr + currLengths[d] - 1)])
        }
        else{
          effLengths[d] = 0
        }
        ctr = ctr + currLengths[d]
      }
    }
    else{
      effLengths = currLengths
    }
    if(length(currLengths) == 0){ return(matrix(0, 0, len+1))}
    ret[,j] = rep(unlist(levels[[j-1]]), effLengths)
    prevLengths = effLengths
  }
  isUnique = apply(ret, MARGIN = 1, function(b) {length(unique(b)) == len + 1})
  ret[isUnique, , drop = FALSE]    
}