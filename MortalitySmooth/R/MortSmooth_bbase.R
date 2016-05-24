MortSmooth_bbase <-
function(x, xl, xr, ndx, deg){
  ## Input:
  ## x   = abcissae of data
  ## xl  = left boundary
  ## xr  = right boundary
  ## ndx = number of internal knots -1
  ##       or number of internal intervals
  ## deg = degree of the splines
  
  ## Output:
  ## B = matrix with the B-spline basis
  
  ## distance between knots
  dx <- (xr - xl) / ndx
  ## One needs (ndx+1) internal knots and 
  ## deg knots on both right and left side
  ## in order to joint all the B-splines
  knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
  ## Truncated deg-th power functions
  ## equally-spaced at given knots along x
  P <- outer(x, knots, MortSmooth_tpower, deg)
  ## number of B-splines which equal to the number of knots
  n <- dim(P)[2]
  ## in the numerator we have the matrix
  ## of deg+1 differences for each knots
  ## this matrix is rescaled by 
  ## (deg! * dx^deg) == (gamma(deg + 1) * dx ^ deg)
  D <- diff(diag(n), diff = deg + 1) /
    (gamma(deg + 1) * dx ^ deg)
  ## the matrix of differences is used to compute B-splines
  ## as differences of truncated power functions
  ## in P %*% t(D)
  ## the last factor is (-1) ^ (deg + 1)
  B <- (-1) ^ (deg + 1) * P %*% t(D)
  B
}
