MortSmooth_BWB <-
function(RTBx, RTBy, nbx, nby, W){
  ## Input:
  ## RTBx: row tensor of the B-spline basis for x
  ## RTBy: row tensor of the B-spline basis for y
  ## nbx: number of B-spline basis for x
  ## nby: number of B-spline basis for y
  ## W: matrix of weights
  
  ## Output:
  ## BWB: a matrix of inner product of a matrix and a sparse weight matrix
  
  BWB <- t(RTBx)%*%W%*%RTBy
  BWB <- array(c(BWB),c(nbx,nbx,nby,nby))
  BWB <- matrix(c(aperm(BWB,c(1,4,2,3))),ncol=nbx*nby)
  return(BWB)
}
