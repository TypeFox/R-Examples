stationary <- function(tau,k,G2,IPI){

  # converto to Pi
  Pi = matrix(0,k,k);
  Pi[IPI] = exp(G2%*%tau)
  Pi = Pi/rowSums(Pi)
  # compute stationary and derivative
  St = Pi;
  D0 = array(0,c(k,k,k^2))
  j = 0;
  for(c2 in 1:k) for(c1 in 1:k){
      j = j+1
      D0[c1,c2,j] = 1
  }
  D = D0
  for(it in 1:10000){
    for(j in 1:k^2){
      D[,,j] = D[,,j]%*%Pi+St%*%D0[,,j]
    }  
    St = St%*%Pi
  }
  la = colMeans(St)
  d1 = NULL
  for(j in 1:(k^2)) d1 = rbind(d1,colMeans(D[,,j]))
  Om = diag(Pi[1,])-Pi[1,]%o%Pi[1,];
  for(c in 2:k){
    piv = Pi[c,]
    Om = blkdiag(Om,diag(piv)-piv%o%piv)
  }
  d0 = t(G2)%*%Om
  d1 = d0%*%d1[IPI,]
  out = list(d0=d0,d1=d1)
  out
}