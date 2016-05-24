linprogPD <- function(x0, A, b, epsilon, pdtol=1e-3, pdmaxiter=50) {
  ## Solves
  ## min ||x||_1   subject to  ||Ax-b||_\infty <= epsilon
  ## Adapted from Matlab code for Dantzig Selector by J. Romberg
  N <- length(x0)
  x0 <- matrix(x0, nrow=N)
  b <- matrix(b, nrow=N)
  
  alpha=0.01;
  beta=0.5;
  mu = 10;

  gradf0 <- matrix(rep(c(0,1), each=N), nrow=2*N)

  if (max(abs(A%*%x0 - b)) > epsilon) {
    stop("Infeasible starting point!")
  }

  x <- x0
  u <- 0.95*abs(x0) + 0.1*max(abs(x0))

  Atr <- A%*%x - b

  fu1 <- x - u
  fu2 <- -x - u
  fe1 <- Atr - epsilon
  fe2 <- -Atr - epsilon
  lamu1 <- -1/fu1
  lamu2 <- -1/fu2
  lame1 <- -1/fe1
  lame2 <- -1/fe2

  AtAv <- t(A)%*%(lame1 - lame2)

  sdg <- - sum(c(fu1, fu2, fe1, fe2)* c(lamu1, lamu2, lame1, lame2))

  tau <- mu*(4*N)/sdg

  rdual <- gradf0 + c(lamu1 - lamu2 + AtAv, -lamu1 - lamu2)
  rcent <- -c(lamu1*fu1, lamu2*fu2, lame1*fe1, lame2*fe2) - 1/tau
  resnorm <- sqrt(sum(rdual^2,rcent^2))

  pditer <- 0
  done <- (sdg < pdtol) | (pditer >= pdmaxiter)

  while(!done) {
    w2 <- -1 - (1/fu1 + 1/fu2)/tau

    sig11 <- -lamu1/fu1 - lamu2/fu2
    sig12 <- lamu1/fu1 - lamu2/fu2
    siga <- -(lame1/fe1 + lame2/fe2)
    sigx <- sig11 - sig12^2/sig11

    w1 <- -( t(A)%*%(1/fe2 - 1/fe1) + 1/fu2 - 1/fu1)/tau
    w1p <- w1 - (sig12/sig11)*w2

    Hp <- t(A)%*%diag(as.vector(siga))%*%A + diag(as.vector(sigx))

    dx <- solve(Hp, w1p)

    if (rcond(Hp) < 1e-14) {
      warning("Ill conditioned matrix.  Previous iterate matrix returned! (May increase perturb/lambda.)")
      xp <- x
      return(xp)
    }
    AtAdx <- A%*%dx

    du <- w2/sig11 - (sig12/sig11)*dx

    dlamu1 <- -(lamu1/fu1)*(dx-du) - lamu1 - 1/(fu1*tau)
    dlamu2 <- -(lamu2/fu2)*(-dx-du) - lamu2 - 1/(fu2*tau)

    dlame1 <- -(lame1/fe1)*(AtAdx) - lame1 - 1/(fe1*tau)
    dlame2 <- (lame2/fe2)*(AtAdx) - lame2 - 1/(fe2*tau)

    AtAdv <- t(A)%*%(dlame1 - dlame2)


    iu1 <- dlamu1 < 0
    iu2 <- dlamu2 < 0
    ie1 <- dlame1 < 0
    ie2 <- dlame2 < 0
    ifu1 <- (dx-du) > 0
    ifu2 <- (-dx-du) > 0
    ife1 <- AtAdx > 0
    ife2 <- AtAdx < 0

    smax <- min( -lamu1[iu1]/dlamu1[iu1], -lamu2[iu2]/dlamu2[iu2], -lame1[ie1]/dlame1[ie1], -lame2[ie2]/dlame2[ie2], -fu1[ifu1]/(dx[ifu1] - du[ifu1]), -fu2[ifu2]/(-dx[ifu2] -du[ifu2]), -fe1[ife1]/AtAdx[ife1], -fe2[ife2]/( - AtAdx[ife2])   )
    smax <- min(1,smax)

    s <- 0.99*smax

    suffdec <- FALSE
    backiter <- 0

    while(!suffdec) {
      xp <- x + s*dx
      up <- u + s*du
      Atrp <- Atr + s*AtAdx
      AtAvp <- AtAv+s*AtAdv

      fu1p <- fu1 + s*(dx - du)
      fu2p <- fu2 + s*(-dx-du)
      fe1p <- fe1 + s*AtAdx
      fe2p <- fe2 + s*(-AtAdx)
      
      lamu1p <- lamu1 + s*dlamu1
      lamu2p <- lamu2 + s*dlamu2
      lame1p <- lame1 + s*dlame1
      lame2p <- lame2 + s*dlame2

      rdp <- gradf0 + c(lamu1p - lamu2p + AtAvp, -lamu1p - lamu2p)
      rcp <- -c(lamu1p*fu1p, lamu2p*fu2p, lame1p*fe1p, lame2p*fe2p) - 1/tau
      suffdec <- sqrt( sum(rdp^2, rcp^2)) < (1 - alpha*s)*resnorm

      s <- beta*s
      backiter <- backiter+1

      if (backiter > 32) {
        warning("Backtracking stuck.  Previous iterate matrix returned!")
        xp <- x
        return(xp)
      }
    }
    x <- xp
    u <- up
    Atr <- Atrp
    AtAv <- AtAvp
    fu1 <- fu1p
    fu2 <- fu2p
    fe1 <- fe1p
    fe2 <- fe2p
    lamu1 <- lamu1p
    lamu2 <- lamu2p
    lame1 <- lame1p
    lame2 <- lame2p

    sdg <- -sum( c(fu1, fu2, fe1, fe2)*c(lamu1, lamu2, lame1, lame2))
    tau <- mu*(4*N)/sdg
    rdual <- rdp
    rcent <- -c(fu1, fu2, fe1, fe2)*c(lamu1, lamu2, lame1, lame2) - 1/tau
    resnorm <- sqrt(sum(rdual^2, rcent^2))

    pditer <- pditer+1
    done <- (sdg < pdtol) | (pditer >= pdmaxiter)
  }

  return(xp)
}
