GmedianCov <- function(X, init=NULL, scores=2, gamma=2, gc=2, alpha=0.75, nstart=1){
  ### Computation of the Geometric covariation matrix 
  ### with averaged stochastic gradient algorithms
  ### input : X   (n x p matrix, n observations, in dimension p)
  ### output : (geometric) median (1 x p numeric vector) and (geometric) median covariation matrix (p x p)  
  ### require library(rARPACK)
  Gmed.est = Gmedian(X,init=init,gamma=gamma,alpha=alpha,nstart=nstart)
  GMCM.est = MedianCovMatRow_rcpp(X,Gmedian=Gmed.est,gamma=gc,alpha=alpha,nstart=nstart)
  if (scores==FALSE){ 
  return(list(median = Gmed.est,covmedian=GMCM.est))
  }
  else {
    ### Computation of the eigenvectors and scores 
    vectors <- RSpectra::eigs_sym(GMCM.est, scores)$vectors
    scores = sweep(X,2,Gmed.est)%*%vectors
    return(list(median=Gmed.est,covmedian=GMCM.est,scores=scores,vectors=vectors))
  }
}

WeiszfeldCov <- function(X, weights=NULL, scores=2, epsilon=1e-08, nitermax = 100){
  ### Computation of the Geometric covariation matrix 
  ### output : (geometric) median (1 x p numeric vector) and (geometric) median covariation matrix (p x p)  
  ### require library(rARPACK)
  X <- as.matrix(X)
  n <- nrow(X)
  if (is.null(weights)) poids <- rep(1,n) 
  Wmed.est <- Weiszfeld_rcpp(X,poids,epsilon=epsilon,nitermax=nitermax)
  WMCM.est <- MedianCovMatW_rcpp(X,Wmed.est$median,poids,epsilon=epsilon,nitermax=nitermax)
  if (scores==FALSE){ 
    return(list(median = Wmed.est$median, covmedian=WMCM.est$median, iterm = Wmed.est$iter, itercov = WMCM.est$iter))
  }
  else {
    ### Computation of the eigenvectors and scores 
    vectors <- RSpectra::eigs_sym(WMCM.est$median, scores)$vectors
    vscores = sweep(X,2,Wmed.est$median)%*%vectors
    return(list(median=Wmed.est$median, covmedian=WMCM.est$median, scores=vscores, vectors=vectors, iterm = Wmed.est$iter, itercov = WMCM.est$iter))
  }
}
