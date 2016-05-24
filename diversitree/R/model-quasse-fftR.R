make.branches.quasse.fftR <- function(control) {
  nx <- control$nx
  dx <- control$dx
  r <- control$r
  dt.max <- control$dt.max
  tc <- control$tc

  f.hi <- make.pde.quasse.fftR(nx*r, dx/r, dt.max, 2L)
  f.lo <- make.pde.quasse.fftR(nx,   dx,   dt.max, 2L)
  combine.branches.quasse(f.hi, f.lo, control)
}

make.pde.quasse.fftR <- function(nx, dx, dt.max, nd) {
  function(y, len, pars, t0) {
    padding <- pars$padding
    ndat <- length(pars$lambda)
    
    nt <- as.integer(ceiling(len / dt.max))
    dt <- len / nt
    if ( !(length(y) %in% (nd * nx)) )
      stop("Wrong size y")
    if ( length(pars$lambda) != length(pars$mu) ||
         length(pars$lambda) > (nx-3) )
      stop("Incorrect length pars")
    if ( pars$diffusion <= 0 )
      stop("Invalid diffusion parameter")

    if ( !is.matrix(y) )
      y <- matrix(y, nx, nd)

    ans <- quasse.integrate.fftR(y, pars$lambda, pars$mu, pars$drift,
                                 pars$diffusion, nt, dt, nx, ndat, dx,
                                 padding[1], padding[2])
    ## Do the log compensation here, to make the careful calcuations
    ## easier later.
    q <- sum(ans[,2]) * dx
    ans[,2] <- ans[,2] / q
    list(log(q), ans)
  }
}

## Note that the sign of drift here needs inverting.  (see
## src/quasse-eqs-fftC.c:qf_setup_kern() for the equivalent flip).
quasse.integrate.fftR <- function(vars, lambda, mu, drift, diffusion,
                                  nstep, dt, nx, ndat, dx, nkl, nkr) {
  kern <- fftR.make.kern(-dt * drift, sqrt(dt * diffusion),
                         nx, dx, nkl, nkr)
  fy <- fft(kern)  
  for ( i in seq_len(nstep) ) {
    vars <- fftR.propagate.t(vars, lambda, mu, dt, ndat)
    vars <- fftR.propagate.x(vars, nx, fy, nkl, nkr)
  }
  vars
}

fftR.make.kern <- function(mean, sd, nx, dx, nkl, nkr) {
  kern <- rep(0, nx)
  xkern <- (-nkl:nkr)*dx
  ikern <- c((nx - nkl + 1):nx, 1:(nkr + 1))
  kern[ikern] <- normalise(dnorm(xkern, mean, sd))
  kern
}

fftR.propagate.t <- function(vars, lambda, mu, dt, ndat) {
  i <- seq_len(ndat)
  r <- lambda - mu
  z <- exp(dt * r)
  e0 <- vars[i,1]
  d0 <- vars[i,-1]
  vars[i,1] <- (mu + z*(e0 - 1)*mu - lambda*e0) /
    (mu + z*(e0 - 1)*lambda - lambda*e0)
  dd <- (z * r * r)/(z * lambda - mu + (1-z)*lambda*e0)^2
  vars[i,-1] <- dd * d0
  vars
}

fftR.propagate.x <- function(vars, nx, fy, nkl, nkr) {
  ifft <- function(x) fft(x, inverse=TRUE)
  f <- function(z, fy, n) {
    nx <- length(z)
    for ( i in seq_len(n) )
      z <- Re(ifft(fft(z) * fy))/nx
    z
  }
  vars.out <- Re(apply(apply(vars, 2, fft) * fy, 2, ifft))/nx
  ndat <- nx - (nkl + 1 + nkr)
  i.prev.l <- 1:nkl
  i.prev.r <- (ndat-nkr+1):ndat
  i.zero <- (ndat+1):nx
  vars.out[c(i.prev.l, i.prev.r),] <- vars[c(i.prev.l, i.prev.r),]
  vars.out[i.zero,] <- 0
  vars.out[vars.out < 0] <- 0
  vars.out
}

