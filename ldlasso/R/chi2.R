chi2 <-
function( Xa, Y ){

  X0 <- Xa[Y==0,]
  X1 <- Xa[Y==1,]
  n0 <- dim(X0)[1]
  n1 <- dim(X1)[1]
  f0 <- colSums(X0)/2/n0
  f1 <- colSums(X1)/2/n1

  OR <- (f1/(1-f1))/(f0/(1-f0))

  O11 <- 2*n0*(1-f0)
  O12 <- 2*n0*f0
  O21 <- 2*n1*(1-f1)
  O22 <- 2*n1*f1

  E11 <- n0/(n0+n1)*(O11 + O21)
  E12 <- n0/(n0+n1)*(O12 + O22)

  E21 <- n1/n0*E11
  E22 <- n1/n0*E12

  X2 <- (O11 - E11)^2/E11 + (O12 - E12)^2/E12 + (O21 - E21)^2/E21 + (O22 - E22)^2/E22

  return(X2)


}

