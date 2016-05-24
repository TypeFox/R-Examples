secder <- function(A,k,l){
  n <- dim(A)[1]
  svec <- matrix(0, nrow = n^2, ncol = n )
  scalesens <- matrix(0, nrow = n^2, ncol = n-1)
  d2 <- matrix(0,nrow = n,ncol = n)
  ev <- eigen(A)
  L <- ev$values

  o <- order(Mod(L))
  W <- ev$vectors
  V <- solve(Conj(W))
  V <- t(Conj(V))

  ### have to re-order V only to match matlab
  o<-order(Mod(L))
  V<-V[,rev(o)]

  for(i in 1:n){
    senmat <- Conj(V[,i])  %*%  t(W[,i])

    svec[,i] <- matrix(senmat,nrow = n^2, ncol = 1)
  }

  s1 <- svec[,1]

  for(m in 2:n){
    scalesens[,m-1] <- svec[,m]/(L[1]-L[m])
  }

   vecsum <- t(apply(scalesens,1,sum))

  for(i in 1:n){
    for(j in 1:n){
      x1 <- (l-1)*n+1 + (i-1)
      x2 <- (j-1)*n+1 + (k-1)
      d2[i,j] <- s1[x1]*vecsum[x2] + s1[x2]*vecsum[x1]
    }
  }

  d2 <- Re(d2)
  d2[A==0] <- 0

  return(d2)
}

