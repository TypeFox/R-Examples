###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2014
#
inputna <-
function(y, S, inputvar=10^4) {
#
################################################################################
# CREATE DESIGN MATRIX X, S, FITTED VALUES AND RESIDUALS
#
  # DIMENSIONS AND NAMES
  y <- as.matrix(y)
  nay <- is.na(y)
  k <- ncol(y)
  nk <- colnames(y)
  if(k>1L && is.null(nk)) nk <- paste("y",seq(k),sep="")
#
  # GENERATE S
  S <- mkS(S,y)
#
  # INPUT
  y[nay] <- 0
  for(i in if(ncol(S)==k) seq(k) else cumsum(c(1,rev(seq(k)[-1]))))
    S[,i][is.na(S[,i])] <- max(S[,i],na.rm=TRUE)*inputvar
  S[is.na(S)] <- 0
#
  cbind(y,S)
}