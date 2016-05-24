mse.saery.boot.AR1.typeIII <-
function(X, D, md, beta, sigma2edi, sigmau1, sigmau2, rho, Fsig){
  
  g1g2g3 <- g123.AR1(X, D, md, sigma2edi, sigmau1, sigmau2, rho, Fsig)
  
  return( g1g2g3[[1]] + g1g2g3[[2]] + 2*g1g2g3[[3]] )
  
}
