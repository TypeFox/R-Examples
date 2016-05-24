## Mixed integer/real optimisation.
do.mle.search.mixed <- function(func, x.init, control, lower, upper) {
  control <- modifyList(list(reltol=.Machine$double.eps^0.25),
                        control)
  is.integer <- control$is.integer

  if ( is.null(is.integer) )
    stop("Must specify which are integers!")
  
  check.bounds(lower, upper, x.init)
  if ( any(is.finite(lower) | is.finite(upper)) )
    func2 <- invert(boxconstrain(func, lower, upper))
  else
    func2 <- invert(func)

  ans <- subplex.mixed(func2, x.init, is.integer,
                       tolerance=control$reltol,
                       scale=control$parscale,
                       verbose=control$verbose)
  ans$value <- -ans$value
  names(ans)[names(ans) == "value"] <- "lnLik"
  ans
}

## Reexpress lik as a "mixed" function....
make.mixed <- function(f, i, xmid, dx) {
  n <- length(i)
  xmid <- check.par.length(xmid, n)
  dx   <- check.par.length(dx, n)

  
  if ( f.is.constrained <- is.constrained(f) ) {
    f.orig <- f
    f <- attr(f, "func")
    i <- match(argnames(f.orig)[i], argnames(f))
  }
  
  ret <- function(x, ..., strict=TRUE) {
    if ( strict )
      check.integer(x[i])
    x[i] <- xmid + dx * x[i]
    f(x, ...)
  }
  class(ret) <- c("mixed", class(f))

  if ( f.is.constrained ) {
    constrain(ret, formulae=attr(f.orig, "formulae"))
  } else {
    ret
  }
}

find.mle.mixed <- function(func, x.init, method, ...) {
  NextMethod("find.mle", method="mixed")
}

trans.mixed <- function(x, i, xmid, dx, to.integer=FALSE) {
  n <- length(i)
  xmid <- check.par.length(xmid, n)
  dx   <- check.par.length(dx, n)

  if ( to.integer )
    x[i] <- round((x[i] - xmid)/dx)
  else
    x[i] <- xmid + dx * x[i]
  x
}

subplex.mixed <- function(func, x.init, is.integer, max.eval=10000,
                          scale=NULL, tolerance=.Machine$double.eps,
                          verbose=FALSE) {
  sign2 <- function(x) (-2*(x <= 0) + 1)
  constrain2 <- function(f, p, i) {
    function(x) {
      p[i] <- x
      f(p)
    }
  }

  psi <- 0.25
  omega <- .1

  nx <- length(x.init) # dimension
  stopifnot(nx >= 2)

  is.integer <- sort(check.integer(is.integer))
  if ( length(is.integer) == 0 || any(is.integer > nx) ||
      any(is.integer < 1) || any(duplicated(is.integer)) )
    stop("Invalid 'is.integer' argument")
  is.real <- setdiff(seq_len(nx), is.integer)

  nx.real <- length(is.real)
  stopifnot(nx.real > 2)
  
  if ( is.null(scale) ) {
    scale <- rep(.1, nx)
    scale[is.integer] <- 1
  } else if ( length(scale) == 1 ) {
    scale <- rep(scale, nx)
    scale[is.integer] <- 1 # Enforce.
  } else if ( length(scale) != nx ) {
    stop("Invalid scale argument")
  }
  scale[is.integer] <- check.integer(scale[is.integer])
  
  subspc.min <- min(2, nx.real)
  subspc.max <- min(5, nx.real)

  cur.x <- x.init
  cur.y <- func(x.init)
  ne <- 1

  delta.x <- step.size <- scale

  flag <- 0
  rep <- 0

  repeat {
    rep <- rep + 1
    prev.x <- cur.x

    delta.x.real <- abs(delta.x[is.real])
    order.x.real <- order(delta.x.real, decreasing=TRUE)

    dims <- partition2(nx.real, order.x.real, delta.x.real,
                       subspc.min, subspc.max)

    ## Simplex loop.
    for ( i in seq_along(dims) ) {
      if ( verbose )
        cat("Subplex...\n")
      j <- is.real[dims[[i]]]
      res <- tryCatch(simplex(constrain2(func, cur.x, j), cur.x[j],
                              step.size[j], max.eval-ne, y.init=cur.y),
                      simplexPrecisionError=function(e) e$ret)
      cur.x[j] <- res$par
      cur.y    <- res$val
      ne <- ne + res$ne
      flag <- res$flag
      
      if ( flag != 0 )
        break # will return
    }

    if ( flag != 0 )
      break # will return

    ## Update integer values, independently.
    for ( j in is.integer ) {
      if ( verbose )      
        cat("Integers...\n")
      res <- minimise.int1d(constrain2(func, cur.x, j), cur.x[j],
                            cur.y, control=list(step=scale[j]))
      cur.x[j] <- res$par
      cur.y    <- res$val
      ne <- ne + res$count
    }
    
    delta.x <- cur.x - prev.x
    delta.x.rel <- pmax(abs(delta.x), abs(step.size)*psi)/
      pmax(abs(cur.x),1)
      
    if ( any(delta.x.rel[is.real] <= tolerance) ) {
      break # will return
    } else {
      ## in this case, we rescale a little, then go back simplex
      if ( length(dims) > 1 )
        step.factor <- min(max(sum(abs(delta.x))/sum(abs(step.size)),
                               omega), 1/omega)
      else
        step.factor <- psi
      step.size <- abs(step.size * step.factor) * sign2(delta.x)
    }
  }

  msg <- list("number of function evaluations exceeds `maxit'",
              NULL,
              "limit of machine precision reached",
              "fstop reached")[[flag+2]]
  list(par=cur.x, value=cur.y, count=ne, convergence=flag,
       message=msg)
}

