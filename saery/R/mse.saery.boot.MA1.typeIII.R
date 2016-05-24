mse.saery.boot.MA1.typeIII <-
function(X, D, md, beta, sigma2edi, sigmau1, sigmau2, theta, Fsig){
  
  g1g2g3 <- g123.MA1(X, D, md, sigma2edi, sigmau1, sigmau2, theta, Fsig)
  
  return( g1g2g3[[1]] + g1g2g3[[2]] + 2*g1g2g3[[3]] )
  
}
