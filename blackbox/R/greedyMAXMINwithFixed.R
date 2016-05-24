greedyMAXMINwithFixed <- function(rownamedarray, finallength, scales, 
                                  fixedNbr=NA ## numeric: number of final rows of rownamedarray, that should be retained in result
                                  ) { # return rownames!
  ## the fixed points should be the last fixedNbr rows in rownamedarray !!
  ## derived from greedy MAXMIN algo, Ravi et al 1991 in Rennen 2008
  #finallength includes fixed points
  if (is.na(fixedNbr)) {
    stop.redef("(!) from greedyMAXMINwithFixed : at least one fixed point required")
    ## peut etre lui faire chercher les deux points les plus distants sans stocker la matrice en memoire
  } else if (fixedNbr<1L) {
    stop("fixedNbr must be a positive integer")
  } ## otherwise greedyMAXMINwithFixed would need to compute the matrix of distances between all points
  if (fixedNbr==finallength) {
    return(rownames(rownamedarray))
  } ## else :
  if (nrow(rownamedarray)<=finallength) {
    if (nrow(rownamedarray)<finallength)
      cat("(!) In greedyMAXMINwithFixed(...): size of input array lower than target final size", "\n")
    return(rownames(rownamedarray))
  } ## else :
  kriglength <- blackbox.getOption("kriglength") ## may be NULL
  if ( ( ! is.null(kriglength)) && (nrow(rownamedarray)-finallength)*finallength > kriglength**2) { ## cf size of final distance matrix
    message.redef("(!) From greedyMAXMINwithFixed: will be slow, and may not have enough memory... ")
  }
  localarray <- sweep(rownamedarray,2L,sqrt(scales),FUN=`/`) # t(t(rownamedarray)/sqrt(scales))
  nnr <- nrow(localarray)
  nnc <- ncol(localarray)
  nnf <- nnr-fixedNbr
  fixed <- (nnf+1):nnr
  subset <- localarray[fixed, , drop=FALSE] ## the fixed points
  candidates <- localarray[1:nnf, , drop=FALSE]
  ##now we compute the distance of all subset points with all nonsubset points
  ##(but not non-subset with non-subset)
  distCS <- as.matrix(proxy::dist(candidates,subset))
  while (nrow(subset)<finallength) {
    farthestidx <- which.max(apply(distCS, 1, min))
    farthestv <- candidates[farthestidx, , drop=FALSE]
    subset <- rbind(farthestv, subset) ## fixed points remain in last rows
    newcandidates <- candidates[-farthestidx, , drop=FALSE]
    distCS <- distCS[-farthestidx, , drop=FALSE]
    newdist <- as.matrix(proxy::dist(newcandidates,farthestv)) ## dists of remaining candidates to newly selected point
    distCS <- cbind(distCS, newdist) ## order of cols doesn't matter
    candidates <- newcandidates
  }
  return(rownames(subset))  ##includes fixed points
}
