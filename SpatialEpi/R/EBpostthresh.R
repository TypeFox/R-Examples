EBpostthresh <-
function(Y, E, alpha, beta, Xrow=NULL, rrthresh){
  if (is.null(Xrow)) Xrow <- matrix(rep(1,length(Y)),nrow=length(Y),ncol=1)
  mu <- as.numeric(exp(Xrow %*% beta))
  thresh <- 1-pgamma(rrthresh,alpha+Y,(alpha+E*mu)/mu)
  return(thresh=thresh)
}
