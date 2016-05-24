commutation<-function(m,n){
  K<-matrix(0,m*n, m*n)
  H<-matrix(0,m,n)
  for(i in 1:m){
    for(j in 1:n){ 
      H[i,j]<-1
      K<-K+kronecker(H,t(H))
      H[i,j]<-0
    }
  }
  K
}
