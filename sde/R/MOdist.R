MOdist <- function(x, M=50, rangeval=range(x, na.rm=TRUE, finite = TRUE)){

 N <- NROW(x)
 nSeries = NCOL(x)
 if(nSeries==1)
  return(0)
 if(N<2)
  return(NA)
 
 if(any(is.na(x))){
  message("Interpolating missing data\n")
  x <- zoo(x)
  x <- na.approx(x)
 }

 base <- create.bspline.basis(rangeval=rangeval, norder=10, nbasis=M)
 
 MOmat <- function(X){
  psi <- eval.basis(as.numeric(X), base)
  P <- matrix(NA, M, M)
  for(i in 1:M)
   for(j in 1:M)
    P[i,j] <- sum(psi[1:(N-1),i]*psi[2:N,j] + psi[1:(N-1),j]*psi[2:N,i])/(2*N)
  return(invisible(P))
 }

  
 V <- vector(nSeries, mode="list")
# Calculate all the Markov operators
 for(i in 1:nSeries)
  V[[i]] <- MOmat(x[,i])

 D <- matrix(0, nSeries, nSeries)
 colnames(D) <- colnames(x)
 rownames(D) <- colnames(x)

# Calculate distance matrix
 for(i in 1:(nSeries-1))
  for(j in (i+1):nSeries){
   D[j,i] <- sum(abs(V[[i]]-V[[j]]))  
  }
  
  as.dist(D)
}

