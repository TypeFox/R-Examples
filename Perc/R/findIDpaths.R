#' find all paths of a certain length for an individual 
#' 
#' \code{findIDpaths} identifies all unique win-loss pathways of order \eqn{(len - 1)} beginning at selected \code{ID}
#' 
#' @param conf a matrix of conf.mat class. An N-by-N conflict matrix whose \code{(i,j)}th element is the number of times i defeated j.
#' @param ID a numeric or character vector of length 1. It specifys the subject at the beginning of each pathway.
#' @param len a positive integer of length 1 greater than 2. the length of the win-loss paths to be identified (\eqn{len = order + 1})
#' @return return all win-loss paths of \code{length(len)} beginning at \code{ID}
#' @seealso \code{\link{as.conflictmat}}, \code{\link{findAllPaths}}, \code{\link{countPaths}}
#' @examples
#' confmatrix <- as.conflictmat(sampleEdgelist)
#' path38891 <- findIDpaths(confmatrix, ID = "Kuai", len = 2)
#' @export


findIDpaths = function(conf, ID, len = 2){
  
  # making sure conf is of conf.mat
  if (!("conf.mat" %in% class(conf))){
    conf = as.conflictmat(conf)
  }
  
  # making sure that len is of correct format
  if (len < 2) stop("len should be no smaller than 2.")
  if (len > 6) stop("len should be no greater than 6.")
  if(len %% as.integer(len) != 0) {
    stop("'len' needs to be an integer.")
  }
  
  if (!(as.character(ID) %in% row.names(conf))) {
    stop("ID not found in the conflict matrix. Making sure the ID and the conflict matrix are correct.")
  }
  
  i <- which(row.names(conf) == as.character(ID))
  
#  if(sum(conf[i,] > 0) == 0) return(matrix(0, 0, len+1))
#  # i = 1                            
#  #len = 5
#  levels = list()
#  levels[[1]] = which(conf[i,] > 0)  ###? levels[[1]]?
#  for(j in 2:len){
#    levels[[j]] = lapply(unlist(levels[[j-1]]), function(k) which(conf[k,] > 0))
#  }
#  ret = matrix(0, length(unlist(levels[[len]])), len+1)
#  ret[,1] = i
#  ret[,len+1] = unlist(levels[[len]])
#  if(len == 2){
#    ret[,2] = rep(unlist(levels[[1]]), sapply(levels[[2]], length))
#  }
#  for(j in len:2){
#    #j = 4
#    currLengths = sapply(levels[[j]], length)
#    if(j < len){
#      effLengths = numeric(length(currLengths))
#      ctr = 1
#      for(d in 1:length(effLengths)){
#        if(currLengths[d] != 0){
#          effLengths[d] = sum(prevLengths[ctr:(ctr + currLengths[d] - 1)])
#        }
#        else{
#          effLengths[d] = 0
#        }
#        ctr = ctr + currLengths[d]
#      }
#    }
#    else{
#      effLengths = currLengths
#    }
#    if(length(currLengths) == 0){ return(matrix(0, 0, len+1))}
#    ret[,j] = rep(unlist(levels[[j-1]]), effLengths)
#    prevLengths = effLengths
#  }
#  isUnique = apply(ret, MARGIN = 1, function(b) {length(unique(b)) == len + 1})
#  ret[isUnique,]
  pathMatrix <- IDpaths(conf, i, len)
  if (nrow(pathMatrix) == 0) {
    message(c("no pathways of length ", len, " found starting at "), ID)
    #return(
    #  list(
    #    pathMatrix, 
    #    paste(
    #      c("no pathways found starting at"),
    #      ID,
    #      sep = " "
    #      )
    #    )
    #  )
  } else {
    pathOutputmatrix <- pathMatrix
    for (j in 1:length(pathMatrix)){
      pathOutputmatrix[j] <- row.names(conf)[pathMatrix[j]]
    }
    return(pathOutputmatrix)
  }
}