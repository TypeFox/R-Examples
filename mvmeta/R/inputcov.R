###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2014
#
inputcov <- 
function(sd, cor=0) {
#
################################################################################
# 
  # IF DATA FRAME, TRANSFORM TO A MATRIX/VECTOR
  if(is.data.frame(sd)) sd <- drop(as.matrix(sd))
  if(is.data.frame(cor)) cor <- drop(as.matrix(cor))
#
  # IF sd IS A VECTOR, INTERPRETED AS STD DEV FOR A SINGLE MATRIX
  if(is.vector(sd)) sd <- t(sd)
  k <- ncol(sd)
  m <- nrow(sd)
  # IF UNIVARIATE, SIMPLY RETURN
  if(k==1L) return(sd^2)
#
  # IF cor IS A VECTOR, DEPENDING ON ITS LENGTH, INTERPRETED AS:
  #   1 CORRELATION, THE SAME FOR ALL THE OUTCOMES FOR ALL THE MATRICES
  #   m CORRELATIONS, DIFFERENT BETWEEN BUT CONSTANT WITHIN MATRICES
  #   k(k-1)/2 CORRELATIONS, INDENTICAL FOR ALL THE MATRICES
  if(is.vector(cor)) {
    cor <- if(length(cor)%in%c(1L,m)) matrix(cor,m,k*(k-1)/2) else
      if(length(cor)==k*(k-1)/2) matrix(cor,m,k*(k-1)/2,byrow=TRUE) else
      stop("Dimensions of 'sd' and 'cor' not consistent")
  # IF cor IS A MATRIX, INTERPRETED AS:
  #   THE k x k CORRELATION MATRIX IF m=1
  #   THE m x k(k-1)/2 MATRIX OF CORRELATIONS    
  } else if(is.matrix(cor)) {
    if(all(dim(cor)==k) && m==1L) cor <- t(cor[lower.tri(cor)]) else
      if(any(dim(cor)!=c(m,k*(k-1)/2))) 
        stop("Dimensions of 'sd' and 'cor' not consistent")
  }
  # CHECK CORRELATIONS AND DIMENSIONS
  if(any(cor^2>1)) stop("correlations must be between -1 and 1")
#
  # INPUT
  nk <- colnames(sd)
  vcov <- t(sapply(seq(m), function(i) {
    R <- diag(k)
    R[lower.tri(R)] <- cor[i,]
    R[upper.tri(R)] <- t(R)[upper.tri(R)]
    D <- diag(sd[i,])
    vechMat(D%*%R%*%D)
  }))
#
  if(m==1L) {
    vcov <- xpndMat(vcov)
    dimnames(vcov) <- list(nk,nk)
  } else colnames(vcov) <- vechMat(outer(nk,nk,paste,sep="."))
#
  vcov
}
