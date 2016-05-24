# This function calculates derivatives (if not provided) and degrees of polynomial on each segment
# Variables dictionary:
# x,y - vectors of equal length containing coordinates of the data points
# fx - vector containing derivatives of function on given points.
# If fx is NA, derivatives are estimated within the program.
# shape - character vector specifying which shape properties should be required.
# Should contain only 'monotonicity' and/or 'curvature'
# smoothness - numeric; continuity class require
# tol - tolerance of the program

sps_prep_u <- function(
  x, y, fx,
  shape, maxdeg,
  smoothness, tol) {

  N <- length(x)

  if(N != length(fx) && length(fx) != 1)
    stop("Lengths of vectors provided are inconsistent")

  stopifnot(tol > 0)

  if (is.unsorted(x))
    stop('x must be in increasing order')

  k <- smoothness

  Delta <- matrix(nrow = N + 1)
  Hi <- matrix(nrow = N + 1)

  i <- 2:N
  Hi[i] <- 1 / (x[i] - x[i-1])
  Hi[1] <- Hi[2]; Hi[N+1] <- Hi[N]
  Delta[i] <- (y[i] - y[i-1]) * Hi[i]
  Delta[1] <- 2 * Delta[2] - Delta[3]
  Delta[N+1] <- 2 * Delta[N] - Delta[N-1]

  n <- 2 * k + 1

  # Find derivatives satisfying (if possible) shape constraints
  if(any(is.na(fx))){
    s = "monotonicity" %in% shape + 2 * "curvature" %in% shape
    if (s == 0) fx <- partial_deriv_none(N, Delta, Hi)
    if (s == 1) fx <- partial_deriv_monotonicity(N, Delta, RHO(n, k))
    if (s == 2) fx <- partial_deriv_curvature(N, 1, Delta, Hi, n, k, tol)
    if (s == 3) fx <- partial_deriv_mon_curv(N, 1, Delta, Hi, RHO(n, k), n, k, tol)
  }

  # Find the degrees needed
  deg <- degree_computation_U(fx,k,Delta,Hi,shape,maxdeg,tol,N)

  list(x=x, y=y, k=k, fx=fx, deg=deg, dim=1)
}

degree_computation_U <- function(fx,k,Delta,Hi,shape,maxdeg,tol,N){

  #Initialize values of NBAR,MBAR

  Nv <- rep(2 * k + 1, (N - 1))

  if(length(shape) > 0){
    for(i in 1:(N-1))
      Nv[i] <- degree_segment(fx[i], fx[i+1], Delta[i+1], Hi[i+1],
                              k, shape, maxdeg, tol)
  }
  if(any(Nv == maxdeg)) warning('maximum degree achieved')
  Nv
}


degree_segment <- function(fx0,fx1,Delta,Hi,k,shape,maxdeg,tol){

  # Functions below are used to test monotonicity and concavity/concavity and correspond to [eq. 2.6]

  G1 <- function(T1,T2,T3,T4,T5) T5 * k * T1 + (T4 - T5 * k) * T2 - T4 * T5 * T3
  G2 <- function(T1,T2,T3,T4,T5) T5 * k * T2 + (T4 - T5 * k) * T1 - T4 * T5 * T3
  G3 <- function(T1,T2,T3,T4) k * T1 + (T4 - k) * T2 - T4 * T3
  G4 <- function(T1,T2,T3,T4) k * T2 + (T4 - k) * T1 - T4 * T3

  rhi <- tol*Hi # Set relative tolerance for this program
  nbar <- 2*k +1


  if ("monotonicity" %in% shape) {

    #Check the kind of data along 0010

    KAU <- ifelse(fx0 >= 0 && Delta >= rhi && fx1 >= 0, 1,
                  ifelse(fx0 <= 0 && Delta <= -rhi && fx1 <= 0, -1,
                         0))

    # Find nbar such that constraints are satisfied (unless you reach nbar = maxdeg)

    if (KAU != 0){
      for (n in nbar:maxdeg){
        rho = RHO(n,k)
        if (KAU == 1 &&
            G1(fx0,fx1,Delta,n,rho) < rhi &&
            G2(fx0,fx1,Delta,n,rho) < rhi ||
            KAU == -1 &&
            G1(fx0,fx1,Delta,n,rho) > -rhi &&
            G2(fx0,fx1,Delta,n,rho > -rhi)) {
          break
        }
      }
      nbar <- n
    }
  }

  if("curvature" %in% shape){

    # Check convexity along 01. K01 == 1 means convex, K01 == 2 means concave
    K01 <- ifelse(Delta - fx0 >= rhi && fx1 - Delta >= rhi, 1,
                  ifelse(Delta - fx0 <= -rhi && fx1 - Delta <= -rhi, 2,
                         0))

    # Find nbar s.t spline is convex(concave) along grid line 0010
    if (K01 != 0){
      for(n in nbar:maxdeg){
        if(K01 == 1 &&
           G3(fx0,fx1,Delta,n) > -rhi &&
           G4(fx0,fx1,Delta,n) < rhi ||
           K01 == 2 &&
           G3(fx0,fx1,Delta,n) < rhi &&
           G4(fx0,fx1,Delta,n) > -rhi) {
          break
        }
      }
      nbar<- n
    }
  }
  nbar
}

# This function evaluates the spline at given point. It makes use of the unique
# properties of basis functions b so that values of function and its derivative,
# whether provided by user or computed in program, is preserved.
