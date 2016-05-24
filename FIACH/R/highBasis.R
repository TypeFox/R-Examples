highBasis <-function(N,freq,tr){
  K<-round(2*(N*tr/(freq + 1)))
  n <- t(0:(N-1))
  C <- matrix(0,N,K)
  C[,1]=1/sqrt(N)
  
  for (i in 2:K){C[,i] = sqrt(2/N)*cos(pi*(2*n+1)*(i-1)/(2*N))}
  basis<-C[,2:ncol(C)]
  return(basis)
}