mse.saery.boot.MA1 <-
function(X, D, md, beta, sigma2edi, sigmau1, sigmau2, theta, Fsig, type, B){
  switch(type,
         I = mse.saery.boot.MA1.typeI(X, D, md, beta, sigma2edi, sigmau1, sigmau2, theta, Fsig, B),
         II = mse.saery.boot.MA1.typeII(X, D, md, beta, sigma2edi, sigmau1, sigmau2, theta, B),
         III = mse.saery.boot.MA1.typeIII(X, D, md, beta, sigma2edi, sigmau1, sigmau2, theta, Fsig)
  )
}
