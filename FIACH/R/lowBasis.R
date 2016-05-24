lowBasis <-function(N,freq,tr){
  
  K<-round(2*(N*tr/(freq + 1)))
  n <- t(0:(N-1))
  C <- matrix(0,N,N)
  C[,1]=1/sqrt(N)
  
  for (i in 2:N){C[,i] = sqrt(2/N)*cos(pi*(2*n+1)*(i-1)/(2*N))}
  basis<-C[,K:ncol(C)]
  return(basis)
}
