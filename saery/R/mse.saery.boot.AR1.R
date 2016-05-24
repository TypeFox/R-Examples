mse.saery.boot.AR1 <-
function(X, D, md, beta, sigma2edi, sigmau1, sigmau2, rho, Fsig, type, B){
  switch(type,
         I = mse.saery.boot.AR1.typeI(X, D, md, beta, sigma2edi, sigmau1, sigmau2, rho, Fsig, B),
         II = mse.saery.boot.AR1.typeII(X, D, md, beta, sigma2edi, sigmau1, sigmau2, rho, B),
         III = mse.saery.boot.AR1.typeIII(X, D, md, beta, sigma2edi, sigmau1, sigmau2, rho, Fsig)
  )
}
