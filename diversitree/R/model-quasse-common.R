## Functions that provide useful lambda and mu functions.
sigmoid.x <- function(x, y0, y1, xmid, r)
  y0 + (y1 - y0)/(1 + exp(r * (xmid - x)))
sigmoid2.x <- function(x, y0, y1, xmid, r) 
  y0 + (y1 - y0)/(1 + exp(4 * r * (xmid - x) / (y1 - y0)))
constant.x <- function(x, c) rep(c, length(x))
noroptimal.x <- function(x, y0, y1, xmid, s2)
  y0 + (y1-y0)*exp(-(x - xmid)^2/(2 * s2))
make.linear.x <- function(x0, x1) {
  if ( is.null(x1) ) {
    function(x, c, m) {
      x1 <- length(x) - x0 + 1
      x[seq_len(x0)]  <- x[x0]
      x[x1:length(x)] <- x[x1]
      ans <- m * x + c
      ans[ans < 0] <- 0
      ans
    }
  } else {
    function(x, c, m) {
      x[x < x0] <- x0
      x[x > x1] <- x1
      ans <- m * x + c
      ans[ans < 0] <- 0
      ans
    }
  }
}
stepf.x <- function(x, y0, y1, xmid)
  ifelse(x < xmid, y0, y1)

normalise <- function(x) x / sum(x)

starting.point.quasse <- function(tree, states, states.sd=NULL) {
  p.bd <- starting.point.bd(tree)

  lik.bm <- make.bm(tree, states, states.sd,
                    control=list(method="pruning", backend="C"))
  c(p.bd, diffusion=as.numeric(stats::coef(find.mle(lik.bm, .1))))
}

load.wisdom <- function(file="wisdom") {
  w <- paste(readLines(file), collapse="\n")
  .Call("r_set_wisdom", w, PACKAGE="diversitree")
}

save.wisdom <- function(file="wisdom") {
  w <- .Call("r_get_wisdom", PACKAGE="diversitree")
  write(w, file)
}

## Checking and sanitisation code:
check.f.quasse <- function(f) {
  args <- names(formals(f))
  if ( args[1] != "x" )
    stop("First argument of speciation/extinction function must be x")
  length(args) - 1
}

check.states.quasse <- function(tree, states, states.sd) {
  states <- check.states(tree, states, as.integer=FALSE)

  if ( length(states.sd) == 1 )
    states.sd <- structure(rep(states.sd, length(states)),
                           names=names(states))
  else
    states.sd <- check.states(tree, states.sd, as.integer=FALSE)
  
  list(states=states, states.sd=states.sd)
}  

check.control.quasse <- function(control, tree, states) {
  tree.length <- max(branching.times(tree))
  xr <- diff(range(states))
  xr.mult <- if ( "xr.mult" %in% names(control) )
    control$xr.mult else 5
  defaults <- list(tc=tree.length/10,
                   dt.max=tree.length/1000,
                   nx=1024,
                   dx=xr * xr.mult / 1024,
                   r=4,
                   xmid=mean(range(states)),
                   w=5,
                   method="fftC",
                   tips.combined=FALSE,
                   flags=FFTW.MEASURE, # fftC only
                   atol=1e-6, # mol only
                   rtol=1e-6, # mol only
                   eps=1e-6,  # perhaps scale with dx?
                   verbose=FALSE)

  nx.changed <- "nx" %in% names(control)
  dx.changed <- "dx" %in% names(control)
  control <- if ( is.null(control) )
    defaults else modifyList(defaults, control)
  if ( dx.changed && !nx.changed )
    control$nx <- 2^ceiling(log2(xr * xr.mult / control$dx))
  else if ( nx.changed && !dx.changed )
    control$dx <- xr * xr.mult / control$nx

  ## Eventually, this will contain "mol"
  method <- match.arg(control$method, c("fftC", "fftR", "mol"))

  if ( control$tips.combined && method != "fftC" )
    stop("'tips.combined' can only be used with method 'fftC'")

  if ( control$tc <= 0 || control$tc >= tree.length )
    stop(sprintf("tc must lie in (0, %2.2f)", tree.length))
  if ( log2(control$nx) %% 1 != 0 )
    stop("nx must be a power of two")
  if ( log2(control$r) %% 1 != 0 )
    stop("r must be a power of two")

  rr <- with(control, xmid + c(-1,1) * dx * nx / 2)
  rmin <- min(c(1, -1) * (mean(range(states)) - rr) / (xr / 2))
  if ( rmin - xr.mult < -1e-5 )
    warning("Range does not look wide enough - be careful!")
  else if ( rmin < 2 )
    stop("Range is not wide enough")

  ## These will be passed through to some C code, so type safety is
  ## important.
  ctrl.int <- c("nx", "flags", "verbose")
  ctrl.num <- c("tc", "dt.max", "r", "xmid", "w", "atol", "rtol")
  control[ctrl.int] <- sapply(control[ctrl.int], as.integer)
  control[ctrl.num] <- sapply(control[ctrl.num], as.numeric)

  control
}

## I use a number of elements of pars.
## pars[[i]]${lambda,mu,drift,diffusion,padding}
## pars$tr
expand.pars.quasse <- function(lambda, mu, args, ext, pars) {
  pars.use <- vector("list", 2)
  for ( i in c(1,2) ) {
    x <- list()
    pars.use[[i]] <-
      list(x=ext$x[[i]], # May screw other things up (was $x[i])
           lambda=do.call(lambda, c(ext$x[i], pars[args$lambda])),
           mu=do.call(mu, c(ext$x[i], pars[args$mu])),
           drift=pars[args$drift],
           diffusion=pars[args$diffusion],
           padding=ext$padding[i,],
           ndat=ext$ndat[i],
           nx=ext$nx[i])
  }
  names(pars.use) <- c("hi", "lo")
  pars.use$tr <- ext$tr
  pars.use
}

make.pars.quasse <- function(cache) {
  args <- cache$args

  function(pars) {
    names(pars) <- NULL # Because of use of do.call, strip names

    drift <- pars[args$drift]
    diffusion <- pars[args$diffusion]

    ext <- quasse.extent(cache$control, drift, diffusion)
    ## This would confirm the translation:
    ##   all.equal(ext$x[[1]][ext$tr], ext$x[[2]])

    ## Parameters, expanded onto the extent:
    pars <- expand.pars.quasse(cache$lambda, cache$mu, args, ext, pars)

    check.pars.quasse(pars$hi$lambda, pars$hi$mu, drift, diffusion)

    pars
  }
}

quasse.extent <- function(control, drift, diffusion) {
  nx <- control$nx
  dx <- control$dx
  dt <- control$dt.max
  xmid <- control$xmid
  r  <- control$r
  w <- control$w

  if ( control$method == "mol" ) {
    ndat <- nx*c(r, 1)
    padding <- NULL
  } else {
    mean <- drift * dt
    sd   <- sqrt(diffusion * dt)

    ## Another option here is to compute all the possible x values and
    ## then just drop the ones that are uninteresting?
    nkl <- max(ceiling(-(mean - w * sd)/dx)) * c(r, 1)
    nkr <- max(ceiling( (mean + w * sd)/dx)) * c(r, 1)
    ndat <- nx*c(r, 1) - (nkl + 1 + nkr)

    padding <- cbind(nkl, nkr)
    storage.mode(padding) <- "integer"
  }

  x0.2 <- xmid - dx*ceiling((ndat[2] - 1)/2)
  x0.1 <- x0.2 - dx*(1 - 1/r)

  ## Concatenate the x values, so that the lambda(x), mu(x)
  ## calculations work for both spaces simultaneously.
  x <- list(seq(x0.1, length.out=ndat[1], by=dx/r),
            seq(x0.2, length.out=ndat[2], by=dx))

  tr <- seq(r, length.out=ndat[2], by=r)

  list(x=x, padding=padding, ndat=ndat, tr=tr, nx=c(nx*r, nx))
}

combine.branches.quasse <- function(f.hi, f.lo, control) {
  nx <- control$nx
  dx <- control$dx
  tc <- control$tc
  r <- control$r
  eps <- log(control$eps)
  dt.max <- control$dt.max

  ## This is hacky version of the log compensation.  It also reduces
  ## the stepsize when bisecting a branch.  It doesn't seem to change
  ## much on 
  careful <- function(f, y, len, pars, t0, dt.max) {
    ans <- f(y, len, pars, t0)
    if ( ans[[1]] > eps ) { # OK
      ans
    } else {
      if ( control$method == "fftC" ||
           control$method == "fftR" )
        dt.max <- dt.max / 2 # Possibly needed
      len2 <- len/2
      ans1 <- Recall(f, y,         len2, pars, t0,        dt.max)
      ans2 <- Recall(f, ans1[[2]], len2, pars, t0 + len2, dt.max)
      ans2[[1]][[1]] <- ans1[[1]][[1]] + ans2[[1]][[1]]
      ans2
    }
  }

  ## Start by normalising the input so that eps up there make
  ## sense...
  function(y, len, pars, t0, idx) {
    if ( t0 < tc ) {
      dx0 <- dx / r
      nx0 <- nx * r
    } else {
      dx0 <- dx
      nx0 <- nx
    }

    ## Here, we also squash all negative numbers.
    if ( any(y < -1e-8) )
      stop("Actual negative D value detected -- calculation failure")
    y[y < 0] <- 0
    y <- matrix(y, nx0, 2)
    q0 <- sum(y[,2]) * dx0
    if ( q0 <= 0 )
      stop("No positive D values")
    y[,2] <- y[,2] / q0
    lq0 <- log(q0)
    
    if ( t0 >= tc ) {
      ans <- careful(f.lo, y, len, pars$lo, t0, dt.max)
    } else if ( t0 + len < tc ) {
      ans <- careful(f.hi, y, len, pars$hi, t0, dt.max)
    } else {
      len.hi <- tc - t0
      ans.hi <- careful(f.hi, y, len.hi, pars$hi, t0, dt.max)

      y.lo <- ans.hi[[2]][pars$tr,]
      lq0 <- lq0 + ans.hi[[1]]
      if ( nrow(y.lo) < nx )
        y.lo <- rbind(y.lo, matrix(0, nx - length(pars$tr), 2))

      ## Fininshing up with the low resolution branch...
      ans <- careful(f.lo, y.lo, len - len.hi, pars$lo, tc, dt.max)
    }

    c(ans[[1]] + lq0, ans[[2]])
  }
}
