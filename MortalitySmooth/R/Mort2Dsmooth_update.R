Mort2Dsmooth_update <-
function(x, y, Z, offset,
                                psi2, wei, Bx, By, nbx, nby,
                                RTBx, RTBy, P, a){
  ## Input:
  ## x: abcissae of data
  ## y: ordinate of data
  ## Z: matrix of count response
  ## offset: a matrix of an a priori known component
  ## psi2: overdispersion parameter
  ## wei: a matrix of weights
  ## Bx: B-splines basis for the x-axis
  ## By: B-splines basis for the y-axis
  ## nbx: number of B-splines for the x-axis
  ## nby: number of B-splines for the y-axis
  ## RTBx: tensors product of B-splines basis for the x-axis
  ## RTBy: tensors product of B-splines basis for the y-axis
  ## P: total penalty matrix for both axes
  ## a: coefficients (in a matrix)
  
  ## Output:
  ## a: updated coefficients (in a matrix)
  
  eta <- MortSmooth_BcoefB(Bx, By, a)
  mu <- exp(offset + eta)
  W <- mu
  z <- eta + (1/mu)* (Z - mu)
  z[which(wei==0)] <- 0
  WW <- wei*W
  BWB <- MortSmooth_BWB(RTBx, RTBy, nbx, nby, WW)
  BWz <- MortSmooth_BcoefB(t(Bx),t(By),(WW*z))
  a0 <- solve(BWB + psi2*P, c(BWz))
  a <- matrix(a0, nrow = nbx)
  a    
}
