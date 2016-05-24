# This function calculates partial derivatives and dergree of polynomial (if not
# provided), according to shape required.

# Variables dictionary:
# X - vector of length np, contains x-coordinates of grid points
# Y - vector of length mp, contains y-coordinates of grid points
# Z - matrix of dimentions (np,mp), contains values of function on grid points
# fx - matrix of dimentions (np,mp), contains values of x partial derivative
# fx[i,j] is the value of x partial derivative on point (X[i],Y[j])
#
# fy - matrix of dimentions (np,mp), contains values of y partial derivative
# fy[i,j] is the value of y partial derivative on point (X[i],Y[j])
#
# Fxy - matrix of dimentions (np,mp) containaing values of the mixed partial
# derivative.  Fxy[i,j] is the value of the mixed partial derivative on point
# (X[i], Y[j])
#
# shape - whether to impose constraints that monotonicity or curvature be
# preserved.  These can also be specified for each dimension separately with
# shape.x and shape.y.
#
# max.deg - maximum degree allowed.  If the function reaches the maximum degree,
# than shape might not be preserved.
#
# smoothness - continuity class required
#
# tol - relative tolerance used by program.  This must be greater than the
# machine precision.

sps_prep_bi = function(
  x, y, z,
  fx, fy, fxy,
  shape.x, shape.y,
  max.deg, smoothness, tol)
{
  N <- length(x)
  M <- length(y)

  if (is.unsorted(x) || is.unsorted(y))
    stop('x and y data must be in increasing order')

  k <- smoothness
  stopifnot(nrow(z) == N)
  stopifnot(ncol(z) == M)

  stopifnot(tol > 0)

  # Check the dimensions of the derivatives (if provided)
  if (all(!is.na(fx))) stopifnot(dim(fx) == c(N, M))
  if (all(!is.na(fy))) {
    stopifnot(dim(fy) == c(N, M))
    fy <- t(fy)	# This is convenient from a mathematical point of view.
  }
  if (all(!is.na(fxy))) stopifnot(dim(fxy) == c(N, M))
  derivs <- list(fx=fx, fy=fy, fxy=fxy)

  # Calculate, Hi, Zi, Delta, Theta.
  # Hi, Zi - are reciprocals of difference between grid points
  # Delta (Theta) - Values of slope of chords along grid lines in x (y)
  # dimension

  Delta <- matrix(ncol=M, nrow=N+1)
  Theta <- matrix(ncol=N, nrow=M+1)
  Zi <- vector(length=M+1)
  Hi <- vector(length=N+1)

  i <- 2:N
  Hi[i] <- 1 / (x[i] - x[i-1])
  Hi[1] <- Hi[2]
  Hi[N+1] <- Hi[N]
  Delta[i,] <- (z[i,] - z[i-1,]) * Hi[i]
  Delta[1,] <- 2 * Delta[2,] - Delta[3,]
  Delta[N+1,] <- 2 * Delta[N,] - Delta[N-1,]

  i <- 2:M
  Zi[i] <- 1 / (y[i] - y[i-1])
  Zi[1] <- Zi[2]
  Zi[M+1] <- Zi[M]
  Theta[i,] <- t((z[,i] - z[,i-1]) * Zi[i])
  Theta[1,] <- 2 * Theta[2,] - Theta[3,]
  Theta[M+1,] <- 2 * Theta[M,] - Theta[M-1,]

  # If some partial derivatives were not given, then calculate appropriate
  # values [eq. 4.2]
  if (any(is.na(derivs)))
    derivs <- calc_derivs(
      N, M, derivs, shape.x, shape.y,
      Delta, Theta, Hi, Zi,
      tol, k, max.deg)

  # If the polynomial degrees were not given, then compute appropriate
  # degrees [p.494]
  deg <- calc_degree(
    derivs, k,
    Delta, Hi, Theta, Zi,
    shape.x, shape.y, max.deg, tol, N, M)

  c(list(x=x, y=y, z=z, k=k, dim=2), derivs, deg)
}



# This function calculates partial derivatives according to shape constraints
# placed.
calc_derivs <- function(N,M,derivs,shape.x,shape.y,
                        Delta,Theta,Hi,Zi,tol,k,max.deg)
{
  encode.shape <- function(shape)
    "monotonicity" %in% shape + 2 * "curvature" %in% shape

  sx <- encode.shape(shape.x)
  sy <- encode.shape(shape.y)

  fxy.given <- all(!is.na(derivs$fxy))
  if (!fxy.given)
    derivs$fxy <- matrix(NA, nrow=N, ncol=M)

  n <- 2*k + 1

  # If the fx derivatives are missing, then calculate appropriate values.
  if (any(is.na(derivs$fx))) {
    if (sx == 0) derivs$fx <- partial_deriv_none(N, Delta, Hi)
    if (sx == 1) derivs$fx <- partial_deriv_monotonicity(N, Delta, RHO(n, k))
    if (sx == 2) derivs$fx <- partial_deriv_curvature(N, M, Delta, Hi, n, k, tol)
    if (sx == 3) derivs$fx <- partial_deriv_mon_curv(N, M, Delta, Hi, RHO(n, k), n, k, tol)
    if (length(shape.y) > 0) {
      mod <- DerMod(
        N, M, Delta, Theta, Zi, tol,
        derivs$fx, sx, sy, 'x', derivs$fxy)
      derivs$fx <- mod$fx
      derivs$fxy <- mod$fxy
    }
  }
  # Similarly for fy...
  if (any(is.na(derivs$fy))) {
    if (sy == 0) derivs$fy <- partial_deriv_none(M, Theta, Zi)
    if (sy == 1) derivs$fy <- partial_deriv_monotonicity(M, Theta, RHO(n, k))
    if (sy == 2) derivs$fy <- partial_deriv_curvature(M, N, Theta, Zi, n, k, tol)
    if (sy == 3) derivs$fy <- partial_deriv_mon_curv(M, N, Theta, Zi, RHO(n, k), n, k, tol)
    if (length(shape.x) > 0) {
      mod <- DerMod(
        M, N, Theta, Delta, Hi, tol,
        derivs$fy, sy, sx, 'y', derivs$fxy)
      derivs$fy <- mod$fx
      derivs$fxy <- mod$fxy
    }
  }
  if (!fxy.given)
    derivs$fxy <- mixed_partial_derivs(N, M, derivs, Hi, Zi)

  derivs
}


calc_degree <- function(derivs, k, Delta, Hi, Theta, Zi,
                        shape.x, shape.y, max.deg, tol, N, M)
{
  fx <- derivs$fx; fy <- derivs$fy; fxy <- derivs$fxy

  NBAR <- matrix(nrow=N-1, ncol=M-1)
  MBAR <- matrix(nrow=N-1, ncol=M-1)

  for (j in 1:(M-1)) for (i in 1:(N-1)) {
    # Start with the lowest possible degrees for fulfilling the
    # smoothness requirements.
    NBAR[[i,j]] <- 2*k+1
    MBAR[[i,j]] <- 2*k+1

    # Do the shape constraints in the x direction require a higher
    # degree?
    if (length(shape.x) != 0) {
      degs <- calc_degree_rectangle(
        fx[i,j], fx[i+1,j], fx[i,j+1], fx[i+1,j+1],
        fy[j,i], fy[j,i+1], fy[j+1,i], fy[j+1,i+1],
        fxy[i,j], fxy[i+1,j], fxy[i,j+1], fxy[i+1,j+1],
        Delta[i+1,j], Delta[i+1,j+1], Hi[i+1], Zi[j+1],
        k, shape.x, max.deg, tol,
        nbar=NBAR[i,j], mbar=MBAR[i,j])
      NBAR[i,j] <- degs$nbar
      MBAR[i,j] <- degs$mbar
    }

    # Do the shape constraints in the y direction require a higher
    # degree?
    if (length(shape.y) != 0
        && (NBAR[i,j] < max.deg || MBAR[i,j] < max.deg)) {
      degs <- calc_degree_rectangle(
        fy[j,i], fy[j+1,i], fy[j,i+1], fy[j+1,i+1],
        fx[i,j], fx[i,j+1], fx[i+1,j], fx[i+1,j+1],
        fxy[i,j], fxy[i,j+1], fxy[i+1,j], fxy[i+1,j+1],
        Theta[j+1,i], Theta[j+1,i+1], Zi[j+1], Hi[i+1],
        k, shape.y, max.deg, tol,
        nbar=MBAR[i,j], mbar=NBAR[i,j])
      MBAR[i,j] <- degs$nbar
      NBAR[i,j] <- degs$mbar
    }
  }
  Nv <- apply(NBAR,1,max)
  Mv <- apply(MBAR,2,max)

  if(max(Nv) == max.deg || max(Mv) == max.deg) warning('maximum degree reached')

  list(deg.x=Nv, deg.y=Mv)
}




calc_degree_rectangle <- function(Fx00,Fx10,Fx01,Fx11,Fy00,Fy01,Fy10,Fy11,
                                  Fxy00,Fxy10,Fxy01,Fxy11,D0010,D0111,
                                  Hi,Zi,K,Shape,max.deg,tol,nbar,mbar){

  # Functions below are used to test monotonicity and concavity/concavity and correspond to [eq. 2.6]

  G1 <- function(T1,T2,T3,T4,T5) T5 * K * T1 + (T4 - T5 * K) * T2 - T4 * T5 * T3
  G2 <- function(T1,T2,T3,T4,T5) T5 * K * T2 + (T4 - T5 * K) * T1 - T4 * T5 * T3
  G3 <- function(T1,T2,T3,T4) K * T1 + (T4 - K) * T2 - T4 * T3
  G4 <- function(T1,T2,T3,T4) K * T2 + (T4 - K) * T1 - T4 * T3

  rhi <- tol*Hi # Set relative tolerance for this program

  if ("monotonicity" %in% Shape) {

    #Check the kind of data along 0010


    KAU <- ifelse(Fx00 >= 0 && D0010 >= rhi && Fx10 >= 0, 1,
                  ifelse(Fx00 <= 0 && D0010 <= -rhi && Fx10 <= 0, -1,
                         0))

    # Find nbar such that constraints are satisfied (unless you reach nbar = max.deg)

    CM1 <- FALSE
    if (KAU != 0){
      for (n in nbar:max.deg){
        rho = RHO(n,K)
        if (KAU == 1 && G1(Fx00,Fx10,D0010,n,rho) < rhi &&
            G2(Fx00,Fx10,D0010,n,rho) < rhi ||
            KAU == -1 && G1(Fx00,Fx10,D0010,n,rho) > -rhi &&
            G2(Fx00,Fx10,D0010,n,rho > -rhi)) {
          CM1 <- TRUE
          break
        }
      }
      nbar <- n
    }

    # Check the kind of data along 0111

    KAU <- ifelse(Fx01 >= 0 && D0111 >= rhi && Fx11 >= 0 , 1,
                  ifelse(Fx01 <= 0 && D0111 <= -rhi && Fx11 <= 0, -1,
                         0))

    # Find nbar such that constraints are satisfied (unless you reach nbar = max.deg)

    CM2 <- FALSE
    if (KAU != 0) {
      for (n in nbar:max.deg) {
        rho=RHO(n,K)
        if (KAU == 1 && G1(Fx01,Fx11,D0111,n,rho) < rhi &&
            G2(Fx01,Fx11,D0111,n,rho) < rhi ||
            KAU == -1 && G1(Fx01,Fx11,D0111,n,rho) > -rhi &&
            G2(Fx01,Fx11,D0111,n,rho) > -rhi) {
          CM2 <- TRUE
          break
        }
      }
      nbar <- n
    }

    #Consider the constraints inside the rectangle R (only if constraints on the edges are satisfied)
    if (CM1 && CM2) {

      K0010 <- ifelse(Fx00 >= rhi && D0010 >= rhi && Fx10 >= rhi , 1,
                      ifelse(Fx00 <= -rhi && D0010 <= -rhi &&  Fx10 <= -rhi, -1,
                             0))
      K0111 <- ifelse(Fx01 >= rhi && D0111 >= rhi && Fx11 >= rhi , 1,
                      ifelse(Fx01 <= -rhi && D0111 <= -rhi &&  Fx11 <= -rhi, -1,
                             0))

      KAU = ifelse(K0010 == K0111, K0010, 0)

      # Increment nbar and mbar unless shape constraints are satisfied
      if (KAU != 0) {
        EX <- FALSE

        Fx00p <- Fx00 + (K / Zi / mbar) * Fxy00
        Fx10p <- Fx10 + (K / Zi / mbar) * Fxy10
        D0010p <- D0010 + (K / Zi / mbar) * (Fy10 - Fy00) * Hi
        Fx01m <- Fx01 - (K / Zi / mbar) * Fxy01
        Fx11m <- Fx11 - (K / Zi / mbar) * Fxy11
        D0111m <- D0111 - (K / Zi / mbar) * (Fy11 - Fy01) * Hi

        while(nbar < max.deg || mbar < max.deg){
          if (EX){

            Fx00p <- Fx00 + (K / Zi / mbar) * Fxy00
            Fx10p <- Fx10 + (K / Zi / mbar) * Fxy10
            D0010p <- D0010 + (K / Zi / mbar) * (Fy10 - Fy00) * Hi
            Fx01m <- Fx01 - (K / Zi / mbar) * Fxy01
            Fx11m <- Fx11 - (K / Zi / mbar) * Fxy11
            D0111m <- D0111 - (K / Zi / mbar) * (Fy11 - Fy01) * Hi
          }

          rho <- RHO(nbar,K)
          G1p <- G1(Fx00p,Fx10p,D0010p,nbar,rho)
          G2p <- G2(Fx00p,Fx10p,D0010p,nbar,rho)
          G1m <- G1(Fx01m,Fx11m,D0111m,nbar,rho)
          G2m <- G2(Fx01m,Fx11m,D0111m,nbar,rho)

          # If shape constraints are satisfied, break. For more explanations
          # look at [eq. 2.6]

          if ((KAU == 1 && Fx00p >- rhi && D0010p >- rhi && Fx10p >- rhi &&
               Fx01m >- rhi && D0111m >- rhi && Fx11m >- rhi && G1p < rhi &&
               G2p < rhi && G1m < rhi && G2m < rhi) || (KAU == -1 &&
                                                        Fx00p < rhi && D0010p < rhi && Fx10p < rhi && Fx01m < rhi &&
                                                        D0111m < rhi && Fx11m < rhi && G1p >- rhi && G2p >- rhi &&
                                                        G1m >- rhi && G2m >- rhi)) break

          # If one of the parameters reached maximum, increase the other.
          # If none of them is at maximum, increase degrees alternatively.

          if (nbar == max.deg){
            mbar <- mbar+1; EX <- TRUE
          } else if (mbar == max.deg) {
            nbar <- nbar+1; EX <- FALSE
          } else {
            EX <- !EX
            if (EX) mbar <- mbar+1
            else nbar <- nbar+1
          }
        }
      }
    }
  }


  if("curvature" %in% Shape){

    # Check convexity along 0010. K0010 == 1 means convex, K0010==2 means concave
    K0010 <- ifelse(D0010 - Fx00 >= rhi && Fx10 - D0010 >= rhi, 1,
                    ifelse(D0010 - Fx00 <= -rhi && Fx10 - D0010 <= -rhi, 2,
                           0))

    # Find nbar s.t spline is convex(concave) along grid line 0010
    CC1 <- FALSE
    if (K0010 != 0){
      for(n in nbar:max.deg){
        if(K0010 == 1 &&
           G3(Fx00,Fx10,D0010,n) > -rhi &&
           G4(Fx00,Fx10,D0010,n) < rhi ||
           K0010 == 2 &&
           G3(Fx00,Fx10,D0010,n) < rhi &&
           G4(Fx00,Fx10,D0010,n) > -rhi) {
          CC1 <- TRUE
          break
        }
      }
      nbar<- n
    }

    # Check convexity along 0111. K0010 == 1 means convex, K0010==2 means concave

    K0111 <- ifelse(D0111 - Fx01 >= rhi && Fx11 - D0111 >= rhi, 1,
                    ifelse(D0111 - Fx01 <= -rhi && Fx11 - D0111 <= -rhi, 2,
                           0))

    #Find nbar s.t spline is convex(concave) along grid line 0111
    CC2 <- FALSE
    if(K0111 != 0){
      for(n in nbar:max.deg){
        if(K0111 == 1 &&
           G3(Fx01,Fx11,D0111,n) > -rhi &&
           G4(Fx01,Fx11,D0111,n) < rhi ||
           K0111 == 2 &&
           G3(Fx01,Fx11,D0111,n) < rhi &&
           G4(Fx01,Fx11,D0111,n) > -rhi) {
          CC2 <- TRUE
          break
        }
      }
      nbar <- n
    }

    #Check the data inside rectangle and adjust nbar,mbar accordingly
    if(CC1 && CC2) {

      KAU <- ifelse(K0010 == K0010, K0010, 0)

      if(KAU != 0) {
        EX <- FALSE

        Fx00p <- Fx00 + (K / Zi / mbar) * Fxy00
        Fx10p <- Fx10 + (K / Zi / mbar) * Fxy10
        D0010p <- D0010 + (K / Zi / mbar) * (Fy10 - Fy00) * Hi
        Fx01m <- Fx01 - (K / Zi / mbar) * Fxy01
        Fx11m <- Fx11 - (K / Zi / mbar) * Fxy11
        D0111m <- D0111 - (K / Zi / mbar) * (Fy11 - Fy01) * Hi
        while(nbar < max.deg || mbar < max.deg){
          if (EX){
            Fx00p <- Fx00 + (K / Zi / mbar) * Fxy00
            Fx10p <- Fx10 + (K / Zi / mbar) * Fxy10
            D0010p <- D0010 + (K / Zi / mbar) * (Fy10 - Fy00) * Hi
            Fx01m <- Fx01 - (K / Zi / mbar) * Fxy01
            Fx11m <- Fx11 - (K / Zi / mbar) * Fxy11
            D0111m <- D0111 - (K / Zi / mbar) * (Fy11 - Fy01) * Hi
          }

          G3p <- G3(Fx00p,Fx10p,D0010p,nbar)
          G4p <- G4(Fx00p,Fx10p,D0010p,nbar)
          G3m <- G3(Fx01m,Fx11m,D0111m,nbar)
          G4m <- G4(Fx01m,Fx11m,D0111m,nbar)

          if (KAU == 1 && D0010p - Fx00p > -rhi &&
              Fx10p - D0010p > -rhi && D0111m - Fx01m > -rhi &&
              Fx11m - D0111m > -rhi && G3p > -rhi && G4p < rhi &&
              G3m > -rhi && G4m < rhi || KAU == 2 && D0010p - Fx00p < rhi &&
              Fx10p - D0010p < rhi && D0111m - Fx01m < rhi &&
              Fx11m - D0111m < rhi && G3p < rhi && G4p > -rhi &&
              G3m < rhi && G4m > -rhi) break


          # If one of the parameters reached maximum, increase the other.
          # If none of them is at maximum, increase degrees alternatively.

          if (nbar == max.deg) {
            EX <- TRUE ; mbar <- mbar+1
          } else if (mbar == max.deg) {
            EX <- FALSE ; nbar <- nbar+1
          } else {
            EX <- !EX
            ifelse(EX, mbar <- mbar+1, nbar <- nbar+1)
          }
        }
      }
    }
  }

  list('nbar'=nbar,'mbar'=mbar)
}

