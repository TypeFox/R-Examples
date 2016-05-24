Mort2Dsmooth_se <-
function(RTBx, RTBy, nbx, nby, BWB.P1){
  ## Input:
  ## RTBx: row tensor of the B-spline basis for x
  ## RTBy: row tensor of the B-spline basis for y
  ## nbx: number of B-spline basis for x
  ## nby: number of B-spline basis for y
  ## W: matrix of weights
  ## BWB.P1: inverse of B'WB + P
  ## Output:
  ## SE: a matrix of standard errors
  ##     for the linear predictor term 
  SE <- array(BWB.P1, c(nbx,nby,nbx,nby))
  SE <- aperm(SE, c(1,3,2,4))
  SE <- matrix(SE, c(nbx^2, nby^2))
  Dim <- c(nrow(RTBx), nrow(RTBy)) 
  SE <- matrix(sqrt(RTBx%*%SE%*%t(RTBy)), Dim)
  SE
}
