updateList <- function(oldlist, newlist){
  ## append newlist to oldlist, deleting members of oldlist with the
  ## same names as members of newlist
  if(length(newlist) == 0)
    return(oldlist)
  dup <- match(names(newlist), names(oldlist), 0)
  dup <- dup[dup > 0]
  if(length(dup) > 0)
    c(oldlist[-dup], newlist)
  else c(oldlist, newlist)
}

updateColumns <- function(oldmat, newmat){
  ## Merge oldmat and newmat, with same-named columns from newmat merged
  ## into the corresponding columns of oldmat.
  if(is.null(newmat)){
    if(is.null(oldmat)) stop("both args NULL")
    newmat <- oldmat
  }
  if(!is.matrix(newmat))
    stop("newmat is not a matrix")
  if(is.null(colnames(newmat)))
    stop("newmat has no column names to merge with")
  
  if(is.null(oldmat)) oldmat <- newmat

  matchvec <- match(colnames(oldmat), colnames(newmat), nomatch = 0)
  if(any(matchvec != 0)){
    dupmat <- mergeSeries(oldmat[, matchvec > 0, drop=F], 
		      newmat[, matchvec, drop=F])
    boink <- cbind(dupmat, newmat, oldmat, union = TRUE)
  }
  else
    boink <- cbind(newmat, oldmat)
  boink <- boink[,!duplicated(colnames(boink)), drop=F]
  boink[, order(colnames(boink)), drop=F]
}
