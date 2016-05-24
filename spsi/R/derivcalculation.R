# This function computes rho, as defined by Costantini [eq. 2.5], which is used
# for calculating the polynomial degrees if monotonicity is required.
#
#                   (n-k-1)                    (n-k-1)         (k-1)
# rho(n) = n/(n-2k)* SUM((n-1)Cv))/((n/(n-2k)* SUM((n-1)Cv))-2*SUM((n-1)Cv)
#                     v=k                       v=k             v=0
RHO <- function(n, k)
{
  S <- sum(choose(n-1, k : (n-k-1)))
  T <- sum(choose(n-1, 0 : (k-1)))
  n / (n - 2*k) * S / (2*k / (n - 2*k) * S - 2*T)
}


G <- function(del1, del2, rho)
{
  ifelse(del1 * del2 <= 0, 0,
         ifelse(abs(del1) < abs(del2),
                rho * del1 * del2 / (del2 + (rho-1) * del1),
                rho * del1 * del2 / (del1 + (rho-1) * del2)))
}

# This function finds the closest point to x that lies inside [lower, upper].
closest.within <- function(x, lower, upper, tol=0.00001)
{
  #print(unlist(list(lower=lower, upper=upper)))
  stopifnot(all(lower <= upper + tol))
  ifelse(x < lower, lower, min(x, upper))

}

# This function computes first order partial derivatives as a weighted average
# of slopes
partial_deriv_none <- function(N, Delta, Hi)
{
  i <- 1:N
  fx <- (Hi[i]*Delta[i,] + Hi[i+1]*Delta[i+1,]) / (Hi[i] + Hi[i+1])
}

# This function computes first order partial derivatives such that constraint
# (monotonicity) placed can be satisfied.
partial_deriv_monotonicity <- function(N, Delta, rho)
{
  i <- 1:N
  fx <- G(Delta[i,], Delta[i+1,], rho)
}


# This function computes first order partial derivatives such that constraint
# (convexity) placed can be satisfied.  The derivatives are computed as a
# sequence of slopes  fx(i,j), i=0,1,..,np , such that (fx(i,j),fx(i+1,j))
# belongs to C(i,j), i=0,1,...,np-1 , where C(i,j) stands for the
# convexity/concavity region (of the n-degree spline along the slice y=y(j)) in
# (x(i),x(i+1)). This sequence exists if and only if ax(1,i) is not greater than
# ax(2,i) , i=0,1,...,np (see the description of the array AX). If, for some
# index  i1 , it has  ax(1,i1) > ax(2,i1) , we say that  x(i1) is a 'break
# point'. It turns out that the one dimensional slice of the spline is co-convex
# everywhere, except for the intervals (x(i1-1),x(i1+1)) which include the break
# points.  There is some arbitrariness in the choice of derivatives, and partial_deriv_curvature
# select the 'closest' values to the corresponding derivatives of Bessel
# interpolation.
partial_deriv_curvature <- function(N,M,Delta,Hi,n,k,tol)
{
  fx <- matrix(nrow=N, ncol=M)

  # These functions define boundaries of concave/convex region
  G3 <- function(T1,T2) (n * T2 - k * T1) / (n - k)
  G4 <- function(T1,T2) (n * T2 - (n - k) * T1) / k


  for(j in 1:M){
    #  Determine the kind (convex or concave) of data in each subinterval
    #  (x(i),x(i+1)) .
    convex <- classify.shape(Delta[,j], Hi, tol)$convex

    # Find the initial set of boundaries for the derivatives
    i <- 1 : (N-1)
    lower <- ifelse(convex, -Inf, Delta[i+1,j])
    upper <- ifelse(convex,  Delta[i+1,j], Inf)

    lower[N] <- -Inf
    upper[N] <- Inf
    I1 <- 1

    while( I1 < N) {

      # Find the boundaries of convexity/concavity regions basing on previous interval
      for (i in (I1+1):N) {
        if(convex[i-1]) {
          P1 <- G3(upper[i-1],Delta[i,j])
          P2 <- G4(lower[i-1],Delta[i,j])
        } else {
          P1 <- G4(upper[i-1],Delta[i,j])
          P2 <- G3(lower[i-1],Delta[i,j])
        }

        if( i == N)
          break
        # If the intervals are disjoint treat this interval as containing inflection point
        if(P1 > (upper[i] + tol * min(Hi[i],Hi[i+1])) ||
           P2 < (lower[i] - tol * min(Hi[i],Hi[i+1])))
          break

        # Otherwise update boundaries for the derivative

        lower[i] <- max(P1,lower[i])
        upper[i] <- min(P2,upper[i])

      }

      # Fix the value of the derivative on next inflection point and subsequently fix earlier derivatives
      I2 <- i # inflection point

      # fxstar -  corresponding  derivative of Bessel interpolation.
      fxstar <- (Hi[I2] * Delta[I2,j] + Hi[I2+1] * Delta[I2+1,j]) / (Hi[I2]+Hi[I2+1])
      fx[I2,j] <- closest.within(fxstar,P1,P2)

      for(i in (I2-1):I1) {

        fxstar <- (Hi[i]*Delta[i,j] + Hi[i+1] * Delta[i+1,j]) / (Hi[i]+Hi[i+1])

        # Ignore bounds at the inflection point
        if (i == I1 && i != 1) {
          fx[i,j] <- fxstar
          break
        }

        if (convex[i]) {
          P1 <- max(lower[i],G4(fx[i+1,j],Delta[i+1,j]))
          P2 <- min(upper[i],G3(fx[i+1,j],Delta[i+1,j]))
        } else {
          P1 <- max(lower[i],G3(fx[i+1,j],Delta[i+1,j]))
          P2 <- min(upper[i],G4(fx[i+1,j],Delta[i+1,j]))
        }
        fx[i,j] <- closest.within(fxstar,P1,P2)
      }

      I1 <- I2
    }
  }
  fx
}

# This function uses the left-slope (and grid width data) at each point i to
# classify the monotonicity and convexity on each interval [i, i+1].
classify.shape <- function(Delta, Hi, tol)
{
  N <- length(Delta) - 2	# number of intervals
  i <- 1:N

  # Classify [i,i+1] as increasing if the right derivative at i is
  # positive (which is equivalent to the data being increasing from i to
  # i+1).
  increasing <- Delta[i+1] >= 0

  # Classify [i,i+1] as convex if the right derivative at i is greater
  # than the left derivative at i.  Note that this is somewhat arbitrary,
  # in that the classification is determined by behaviour at near i (not
  # i+1).
  convex <- Delta[i+1] >= Delta[i]

  # Some amendments:
  for (i in 1:(N-1)) {
    # If the slope is within the tolerance range of being constant
    # over two adjacent intervals, then treat the function as being
    # monotone and of the same curvature over these intervals.
    if (abs(Delta[i] - Delta[i+1]) <= tol * min(Hi[i], Hi[i+1])) {
      increasing[i] <- increasing[i+1]
      convex[i] <- convex[i+1]
    }

    # If there's a local maximum at i+1, then ensure concavity on
    # both sides of i+1.  (And similarly for local minima.)
    if (increasing[i] != increasing[i+1])
      convex[i] <- convex[i+1]
  }

  list(increasing=increasing, convex=convex)
}


# This function selects first order partial derivatives in a way that preserves
# monotonicity and concavity.
partial_deriv_mon_curv <- function(N, M, Delta, Hi, rho, n, k, tol) {

  fx <- matrix(nrow=N, ncol=M)

  G3 <- function(T1, T2) (n*T2 - k*T1) / (n-k)
  G4 <- function(T1, T2) (n*T2 - (n-k)*T1) / k

  for(j in 1:M) {

    class <- classify.shape(Delta[,j], Hi, tol)
    increasing <- class$increasing
    convex <- class$convex

    i <- 1:(N-1)
    lower <- ifelse(convex == increasing, pmin(0, Delta[i+1,j]), pmin(Delta[i+1,j], n*Delta[i+1,j]/k))
    upper <- ifelse(convex == increasing, pmax(0, Delta[i+1,j]), pmax(Delta[i+1,j], n*Delta[i+1,j] / k))
    lower[N] <- -Inf
    upper[N] <- Inf

    I1 <- 1
    while (I1 < N) {

      # Find the next inflection point, and calculate bounds
      # on derivatives.
      for (i in (I1+1):N) {
        if (increasing[i-1] && convex[i-1]) {
          P1 <- G3(upper[i-1], Delta[i,j])
          P2 <- G4(lower[i-1], Delta[i,j])
        } else if (increasing[i-1] && !convex[i-1]) {
          P1 <- max(0, G4(upper[i-1], Delta[i,j]))
          P2 <- G3(lower[i-1], Delta[i,j])
        } else if (!increasing[i-1] && convex[i-1]) {
          P1 <- G3(upper[i-1], Delta[i,j])
          P2 <- min(0, G4(lower[i-1], Delta[i,j]))
        } else { # if (!increasing[i-1] && !convex[i-1])
          P1 <- G4(upper[i-1], Delta[i,j])
          P2 <- G3(lower[i-1], Delta[i,j])
        }

        if (i == N)
          break

        if (P1 > upper[i] + tol * min(Hi[i], Hi[i+1])
            || P2 < lower[i] - tol * min(Hi[i], Hi[i+1]))
          break

        lower[i] <- max(P1, lower[i])
        upper[i] <- min(P2, upper[i])
      }

      I2 <- i		# inflection point
      fxstar <- G(Delta[I2,j], Delta[I2+1,j], rho)
      fx[I2,j] <- closest.within(fxstar, P1, P2)

      # Select derivatives in a globally consistent way for
      # this group of intervals that contains no interior
      # inflection points.
      for (i in (I2-1):I1) {

        fxstar <- G(Delta[i,j], Delta[i+1,j], rho)

        # Ignore bounds at the inflection point
        if (i == I1 && i != 1) {
          fx[i,j] <- fxstar
          break
        }
        if (convex[i]) {
          P1 <- max(lower[i], G4(fx[i+1,j], Delta[i+1,j]))
          P2 <- min(upper[i], G3(fx[i+1,j], Delta[i+1,j]))
        } else {
          P1 <- max(lower[i], G3(fx[i+1,j], Delta[i+1,j]))
          P2 <- min(upper[i], G4(fx[i+1,j], Delta[i+1,j]))
        }
        fx[i,j] <- closest.within(fxstar, P1, P2)
      }

      I1 <- I2
    }
  }
  fx
}


der_mod_line = function(np,Delta,tol,i,J1,J2,Fx,Shape.x,var,Fxy){
  if(Shape.x == 0){

    Fx[i,J1:J2] <- ( min(Fx[i,J1:J2]) + max(Fx[i,J1:J2]) ) / 2

  }  else if (Shape.x == 1) {
    fxmin <- min(Fx[i,J1:J2])
    fxmax <- max(Fx[i,J1:J2])

    if(fxmin <= 0 && fxmax >= 0) fxstar <- 0

    else {
      mina <- min(abs(fxmin),abs(fxmax))
      maxa <- max(abs(fxmin),abs(fxmax))
      gamma <- sqrt(sqrt(max(0,mina-tol))/maxa)

      fxstar <- ifelse(fxmax > 0, (1 - gamma) * mina + gamma * maxa,
                       (gamma - 1) * mina - gamma * maxa)
    }

    Fx[i,J1:J2] <- fxstar

  } else if(Shape.x == 2) {

    # Alfay (Betay) minimum (maximum) values of derivative on each y line of the grid

    Alfay <- apply(Delta[i:(i+1),J1:J2],2,min)
    Betay <- apply(Delta[i:(i+1),J1:J2],2,max)

    fxstar <- (max(Alfay) + min(Betay)) / 2

    Fx[i,J1:J2] <- closest.within(fxstar,Alfay,Betay)

  } else {

    # Alfay (Betay) minimum (maximum) values of derivative on each y line of the grid

    Alfay <- rep(Inf,J2)
    Betay <- rep(-Inf,J2)

    for(j in J1:J2){
      if((Delta[i,j] * Delta[i+1,j]) <= 0) {
        Alfay[j] <- 0
        Betay[j] <- 0

      } else {
        mina <- min(abs(Delta[i,j]),abs(Delta[i+1,j]))
        maxa <- max(abs(Delta[i,j]),abs(Delta[i+1,j]))
        gamma <- sqrt(sqrt(max(0,mina-tol)/maxa))

        Alfay[j] <- ifelse(Delta[i+1,j] > 0, mina,(gamma - 1)* mina - gamma * maxa)
        Betay[j] <- ifelse(Delta[i+1,j] > 0, (1 - gamma) * mina + gamma * maxa, -mina)

      }
    }
    Fx[i,J1:J2] <- closest.within((max(Alfay[J1:J2]) + min(Betay[J1:J2]))/2,Alfay[J1:J2],Betay[J1:J2])
  }

  if (any(is.na(Fxy))) {
    if(var == 'x') Fxy[i,J1:J2] <- 0
    else Fxy[J1:J2,i] <- 0
  }

  list(Fxy=Fxy, Fx=Fx)
}





# If partial derivatives were not provided and particular shape is required,
# this function modifies derivatives computed earlier to facilitate shape
# preserving interpolation
DerMod = function(np,mp,Delta,Theta,Zi,tol,Fx,Shape.x,Shape.y,var,Fxy){

  if(Shape.y == 1) {
    for(i in 1 : np) {
      J1 <- 1

      while (J1 < mp) {
        if(abs(Theta[J1+1,i]) < (tol * Zi[J1+1])) {

          for(J2 in (J1+1):mp)
            if(abs ( Theta[(J1 + 2),i]) > tol * Zi[(J1 + 2)])
              break

          derivs <- der_mod_line(np,Delta,tol,i,J1,J2,Fx,Shape.x,var,Fxy)
          Fx  <- derivs$Fx
          Fxy <- derivs$Fxy
          J1 <- J2
        } else J1 <- J1+1
      }
    }
  } else if(Shape.y == 2) {
    for(i in 1 : np) {
      J1 <- 1
      while (J1 <= (mp - 2)) {

        for(J2 in (J1+1) : mp)
          if (abs(Theta[J1+1,i] - Theta[(J2+1),i]) >= tol * min(Zi[J1+1],Zi[J2+1]))
            break

        if (J2 >= J1 + 2) {
          derivs <- der_mod_line(np,Delta,tol,i,J1,J2,Fx,Shape.x,var,Fxy)
          Fx = derivs$Fx
          Fxy = derivs$Fxy
        }
        J1 <- J2
      }
    }
  } else {
    for(i in 1:np){
      J1 <- 1
      while(J1 <= mp-1) {

        for(J2 in (J1+1) : mp)
          if (abs(Theta[J1+1,i] - Theta[(J2+1),i]) >= tol * min(Zi[J1+1],Zi[J2+1]))
            break

        if(J2 >= (J1+2) || abs(Theta[J1+1,i]) < tol * Zi[J1+1] ) {
          derivs <- der_mod_line(np,Delta,tol,i,J1,J2,Fx,Shape.x,var,Fxy)
          Fx = derivs$Fx
          Fxy = derivs$Fxy
        }
        J1 <- J2
      }
    }
  }
  list(fx=Fx, fxy=Fxy)
}


# This function computes mixed partial derivatives as quantities based on first
# order partial derivatives. [p. 495]
mixed_partial_derivs <- function(np,mp,derivs,Hi,Zi)
{
  fx <- derivs$fx; fy <- derivs$fy; fxy <- derivs$fxy
  for(j in 2:(mp-1)) {
    for(i in 2:(np-1)) {
      if (!is.na(fxy[i, j]))
        next
      fxy[i,j] <- (Hi[i]^2 * (fy[j,i] - fy[j,i-1])
                   + Hi[i+1]^2 * (fy[j,i+1] - fy[j,i])
                   + Zi[j]^2 * (fx[i,j+1] - fx[i,j])
                   + Zi[j+1]^2 * (fx[i,j+1] - fx[i,j])
      ) / (Hi[i] + Hi[i+1] + Zi[j] + Zi[j+1])
    }
  }

  fxy[which(is.na(fxy))] <- 0
  fxy
}
