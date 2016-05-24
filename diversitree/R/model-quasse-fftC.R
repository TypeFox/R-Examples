FFTW.ESTIMATE <- -1
FFTW.MEASURE <- 0
FFTW.PATIENT <- 1
FFTW.EXHAUSTIVE <- 2

## 8: branches
## TODO: I am not sure that this will pick up on the edge case of tc
## being *exactly* t0 or t0 + len (and what happens if rounding error
## puts this around a break).
## This might be something to check during the checking phase that tc
## does not come within 1e-8 of a node?
make.branches.quasse.fftC <- function(control) {
  check.fftC()
  
  nx <- control$nx
  dx <- control$dx
  r <- control$r
  dt.max <- control$dt.max
  tc <- control$tc
  flags <- control$flags
  
  f.hi <- make.pde.quasse.fftC(nx*r, dx/r, dt.max, 2L, flags)
  f.lo <- make.pde.quasse.fftC(nx,   dx,   dt.max, 2L, flags)
  combine.branches.quasse(f.hi, f.lo, control)
}

make.pde.quasse.fftC <- function(nx, dx, dt.max, nd, flags) {
  ptr <- .Call("r_make_quasse_fft", as.integer(nx), as.numeric(dx),
               as.integer(nd), as.integer(flags),
               PACKAGE="diversitree")
  function(y, len, pars, t0, dt=dt.max) {
    nt <- as.integer(ceiling(len / dt.max))
    dt <- len / nt
    if ( !(length(y) %in% (nd * nx)) )
      stop("Wrong size y")
    if ( any(is.na(y)) )
      stop("Illegal NA values in input")
    if ( length(pars$lambda) != length(pars$mu) ||
         length(pars$lambda) > (nx-3) )
      stop("Incorrect length pars")
    if ( pars$diffusion <= 0 )
      stop("Invalid diffusion parameter")
    
    ans <- .Call("r_do_integrate",
                 ptr, y, pars$lambda, pars$mu,
                 pars$drift, pars$diffusion,
                 nt, dt, pars$padding, PACKAGE="diversitree")

    ## Do the log compensation here, to make the careful calcuations
    ## easier later.
    if ( ncol(ans) > 1 ) {
      q <- sum(ans[,2]) * dx
      ans[,2] <- ans[,2] / q
      list(log(q), ans)
    } else {
      list(0, ans)
    }
  }
}

make.tips.quasse.fftC <- function(control, t, tips) {
  nx <- control$nx
  dx <- control$dx
  r <- control$r
  dt.max <- control$dt.max
  tc <- control$tc
  flags <- control$flags
  
  i <- order(t)
  n <- length(t)
  target <- tips[i]

  tt <- t[i]
  nd <- (n+1):2

  j <- which(tt > tc)[1]
  if ( is.na(j) )
    ## TODO: Fix this bug at some point (not a huge priority, really)
    stop("Please reduce tc (but this is a bug)!")

  nd.hi <- nd[1:j]
  tt.hi <- tt[1:j]
  tt.hi[j] <- tc
  len.hi <- diff(c(0, tt.hi))
  nt.hi <- as.integer(ceiling(len.hi / dt.max))
  dt.hi <- len.hi / nt.hi
  
  nd.lo <- nd[j:n]
  tt.lo <- tt[j:n]
  len.lo <- diff(c(tc, tt.lo))
  nt.lo <- as.integer(ceiling(len.lo / dt.max))
  dt.lo <- len.lo / nt.lo
  
  ptr.hi <- .Call("r_make_quasse_fft", as.integer(nx*r), as.numeric(dx/r),
                  as.integer(nd.hi), as.integer(flags),
                  PACKAGE="diversitree")
  ptr.lo <- .Call("r_make_quasse_fft", as.integer(nx), as.numeric(dx),
                  as.integer(nd.lo), as.integer(flags),
                  PACKAGE="diversitree")

  dx.v <- rep(c(dx/r, dx), c(length(nd.hi)-1, length(nd.lo)))
  
  function(y, pars) {
    if ( is.matrix(y) ) {
      if ( nrow(y) != (nx*r) || ncol(y) != nd.hi[1] )
        stop("Invalid dimension y")
    } else {
      if ( length(y) != (nd.hi[1] * nx * r) )
        stop("Invalid size y")
      y <- matrix(y, nx*r, nd.hi[1])
    }

    if ( pars$hi$diffusion <= 0 || pars$lo$diffusion <= 0 )
      stop("Invalid diffusion parameter")

    ## The columns need to be reversed...
    y.hi <- y[,c(1,rev(2:ncol(y)))]
    
    ret.hi <- .Call("r_do_tips",
                    ptr.hi, y.hi, pars$hi$lambda, pars$hi$mu,
                    pars$hi$drift, pars$hi$diffusion, nt.hi, dt.hi,
                    as.integer(pars$hi$padding),
                    PACKAGE="diversitree")

    y.lo <- rbind(ret.hi[[j]][pars$tr,],
                  matrix(0, nx - length(pars$tr), nd.lo[1]))

    ret.lo <- .Call("r_do_tips",
                    ptr.lo, y.lo, pars$lo$lambda, pars$lo$mu,
                    pars$lo$drift, pars$lo$diffusion, nt.lo, dt.lo,
                    as.integer(pars$lo$padding),
                    PACKAGE="diversitree")

    base <- c(ret.hi[-j], ret.lo)
    lq <- numeric(length(base))
    for ( i in seq_along(base) ) {
      q <- sum(base[[i]][,2]) * dx.v[i]
      base[[i]][,2] <- base[[i]][,2] / q
      lq[i] <- log(q)
    }

    list(target=target, base=base, lq=lq)
  }
}

make.branches.aux.quasse.fftC <- function(control, sampling.f) {
  nx <- control$nx
  dx <- control$dx
  r <- control$r
  dt.max <- control$dt.max
  tc <- control$tc
  flags <- control$flags
  
  pde.hi <- make.pde.quasse.fftC(nx*r, dx/r, dt.max, 1L, flags)
  pde.lo <- make.pde.quasse.fftC(nx,   dx,   dt.max, 1L, flags)
  nx.hi <- nx*r
  nx.lo <- nx

  n <- length(sampling.f)
  e0 <- 1 - unlist(sampling.f)

  function(i, len, pars, idx) {
    if ( i > n )
      stop("No such partition")
    if ( length(len) > 1 )
      stop("Can't do this any more.")

    ndat.hi <- length(pars$hi$lambda)
    npad.hi <- nx.hi - ndat.hi
    npad.lo <- nx.lo - length(pars$lo$lambda)

    y <- matrix(rep(c(e0[i], 0), c(ndat.hi, npad.hi)), nx.hi, 1)
    if ( len < tc ) {
      ans <- pde.hi(y, len,    pars$hi, 0 )[[2]][,1]
    } else {
      y   <- pde.hi(y, tc,     pars$hi, 0 )[[2]][,1]
      y   <- c(y[pars$tr], rep(0, npad.lo))
      ans <- pde.lo(y, len-tc, pars$lo, tc)[[2]][,1]
    }
    ans
  }
}

check.fftC <- function(error=TRUE) {
  ok <- is.loaded("r_make_quasse_fft", "diversitree")
  if ( error && !ok )
    stop("diversitree was built without FFTW support.  ",
         "See Details in ?make.quasse")
  ok
}
