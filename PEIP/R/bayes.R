bayes <-
function( G,mprior,covm,d,covd)
{

  
  covmp=pracma::inv( t(G) %*% pracma::inv(covd) %*% G + pracma::inv(covm)  )
### 
###  This takes care of any lack of symmetry in covmp.
### 
  covmp=(covmp+t(covmp) )/2
  
  covd12=pracma::sqrtm( pracma::inv(covd) )$B
  
  covm12=pracma::sqrtm( pracma::inv(covm) )$B
  
  A=rbind(covd12%*%G,  covm12)

  rhs=c(covd12%*%d, covm12 %*% mprior)
  
  mmap=Ainv(A, rhs, tol=1e-12)
  


  return(list(covmp=covmp,mmap=mmap))

}
