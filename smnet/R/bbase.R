
# Truncated p-th power function
tpower <- function(x, t, p)  (x - t) ^ p * (x > t)

# Construct B-spline basis, not sparse
bbase 	<- 	function(x, xl = min(x), xr = max(x), nseg = 10, deg = 3){
  dx 	<- 	(xr - xl) / nseg
  knot<- 	seq(xl - deg * dx, xr + deg * dx, by = dx)
  P 	<- 	outer(x, knot, tpower, deg)
  n 	<- 	dim(P)[2]
  D 	<- 	diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx ^ deg)
  B 	<- 	(-1) ^ (deg + 1) * P %*% t(D)
  B}


# Construct sparse B-spine basis
b_spline_basis  	<-	function(x, xl = min(x), xr = max(x), nseg = 10, deg = 3){
  dx <- (xr - xl) /nseg
  knots<- seq(xl - deg * dx, xr + deg * dx, by = dx)
  sparseBasis<-splineDesign(knots = knots, x = x, ord = (deg + 1), outer.ok=T, sparse = T)
  make_spam(sparseBasis)
}


cyclic_b_spline_basis<-function(x, xl = min(x), xr = max(x), nseg = 10, deg = 3)
{
  
  dx <- (xr - xl) /nseg
  knots<- seq(xl - deg * dx, xr + deg * dx, by = dx)
  nk <- length(knots)
  k1 <- knots[1]
  xc <- knots[nk - deg + 1]
  knots <- c(k1 - (knots[nk] - knots[(nk - deg  + 1):(nk - 1)]), 
             knots)
  ind <- x > xc
  X1 <- splineDesign(knots = knots, x = x, ord = (deg + 1), outer.ok = TRUE, sparse = T)
  x[ind] <- x[ind] - max(knots) + k1
  if (sum(ind)) {
    X2 <- splineDesign(knots = knots, x = x[ind], ord = (deg + 1), outer.ok = TRUE, sparse = T)
    X1[ind, ] <- X1[ind, ] + X2
  }
  make_spam(X1)
}
