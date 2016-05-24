###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012-2014
#
mkS <- 
function(S, y, narm=NULL, subset=NULL) {
#
################################################################################
# TRANSFORM S IN A MATRIX OF VECTORIZED VCOV MATRICES
#
  # DIMENSIONS
  k <- dim(y)[2]
  m <- dim(y)[1]
#
  # MESSAGE
  mes <- "incorrect dimensions for 'S'"
#
  # IF DATAFRAME
  if(is.data.frame(S)) S <- as.matrix(S)
#
  # IF NUMERIC
  if(is.numeric(S)) {
    # IF JUST A VECTOR, TRANSFORM IN A MATRIX
    if(is.null(dim(S))) S <- as.matrix(S)
    # IF AN ARRAY, TRANSFORM IN A MATRIX
    if(length(dim(S))==3L) S <- t(apply(S,3,vechMat))
    # FINALLY, IF A MATRIX, CHECK DIMENSIONALITY
    if(!is.null(subset)) S <- S[subset,,drop=FALSE]
    if(!is.null(narm)) S <- S[-narm,,drop=FALSE]
    if(dim(S)[1]!=m || !dim(S)[2] %in% c(k,k*(k+1)/2)) stop(mes)    
  }
#
  # IF A LIST
  if(is.list(S)) {
    S <- lapply(S,as.matrix)
    if(!is.null(subset)) S <- S[subset]
    if(!is.null(narm)) S <- S[-narm]
    if(length(S)!=m) stop(mes)
    if(any(sapply(S,dim)!=k)) stop(mes)
    S <- if(k==1L) as.matrix(sapply(S,vechMat)) else t(sapply(S,vechMat))
  }
#
  # NAMES
  rownames(S) <- rownames(y)
  nk <- colnames(y)
  colnames(S) <- if(!is.null(nk)) if(dim(S)[2]==k) nk else
    vechMat(outer(nk,nk,paste,sep=".")) else NULL
#
  S
}

#
