#' Identifies all paths between all pairs of less than or 
#' equal to a certain length
#' 
#' \code{findAllPaths} Identifies all paths length less than or equal 
#' to \code{maxLength} between all pairs of competitors 
#' 
#' @param conf a matrix of conf.mat class. An N-by-N conflict matrix whose \code{(i,j)}th element is the number of times i defeated j.
#' @param maxLength a positive numeric integer indicating the maximum length of paths to identify
#' @return A list of two elements. 
#'  
#'  \item{direct pathways}{direct pathways found in original matrix}
#'  
#'  \item{indirect pathways}{a list of all paths from length 2 to the given length}
#' 
#' @seealso \code{\link{countPaths}} \code{\link{findIDpaths}} \code{\link{transitivity}}
#' @examples
#' # convert an edgelist to conflict matrix
#' confmatrix <- as.conflictmat(sampleEdgelist)
#' # find all paths of legnth 3
#' allp.3 <- findAllPaths(confmatrix, 3)
#' @export


findAllPaths = function(conf, maxLength = 2){
  
  # making sure conf is of conf.mat
  if (!("conf.mat" %in% class(conf))){
    conf = as.conflictmat(conf)
  }
  
  # check if maxLength is correct
  if (maxLength < 2) stop("'maxLength' should be no smaller than 2.")
  if (maxLength > 6) stop("'maxLength' should be no greater than 6.")
  if(maxLength %% as.integer(maxLength) != 0) {
    stop("'maxLength' needs to be an integer.")
  }
  
  allPathsOutput <- allPaths(conf, maxLength)
#  paths = lapply(2:maxLength, FUN = function(l, conf){
#    do.call(rbind, lapply(1:nrow(conf), FUN = IDpaths, conf = conf, l = l))
#  }, conf = conf)
  paths <- allPathsOutput[[2]]
  pathOutput <- paths
  for (i in 1:length(paths)){
    for (j in 1:length(paths[[i]]))
    pathOutput[[i]][j] <- row.names(conf)[paths[[i]][j]]
  }
  pathOutputAll <- list(which(conf > 0, arr.ind = TRUE), pathOutput)
  names(pathOutputAll)[1] <- "direct pathways"
  names(pathOutputAll)[2] <- "indirect pathways"
  return(pathOutputAll)
}

# to do: 
#     -1: #!
#     -2: # name each element in the list with "pathways of length n"