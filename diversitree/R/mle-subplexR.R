## R implementation of subplex.  This exists only for reference, and
## is not exported.  Some of the auxilliary functions (partition,
## simplex) are used by a real/integer mixed optimisation algorithm
## (mle-mixed.R)
do.mle.search.subplexR <- function(func, x.init, control, lower,
                                   upper) {
  control <- modifyList(list(reltol=.Machine$double.eps^0.25,
                             parscale=rep(.1, length(x.init))),
                        control)

  check.bounds(lower, upper, x.init)
  if ( any(is.finite(lower) | is.finite(upper)) )
    func2 <- invert(boxconstrain(func, lower, upper))
  else
    func2 <- invert(func)
  ans <- subplexR(x.init, func2, control)
  ans$value <- -ans$value
  names(ans)[names(ans) == "value"] <- "lnLik"

  if ( ans$convergence != 0 )
    warning("Convergence problems in find.mle (subplex): ",
            tolower(ans$message))
  
  ans
}

subplexR <- function(par, fn, control=list(), hessian=FALSE, ...) {
  sign2 <- function(x) (-2*(x <= 0) + 1)
  constrain2 <- function(f, p, i) {
    function(x) {
      p[i] <- x
      f(p)
    }
  }

  if ( hessian )
    .NotYetUsed("hessian")

  nx <- length(par) # dimension
  if ( nx < 2 )
    stop("Need at least two dimensions for subplexR")

  control <- modifyList(list(maxit=10000, parscale=rep(1, nx),
                             reltol=.Machine$double.eps), control)
  max.eval <- control$maxit
  scale <- control$parscale
  tolerance <- control$reltol
  func <- function(x) fn(x, ...)

  psi <- 0.25
  omega <- .1

  if ( length(scale) == 1 )
    scale <- rep(scale, nx)
  else if ( length(scale) != nx )
    stop("Invalid scale argument")
  
  subspc.min <- min(2, nx)
  subspc.max <- min(5, nx)

  cur.x <- par
  cur.y <- func(par)
  ne <- 1

  delta.x <- step.size <- scale

  flag <- 0
  rep <- 0

  repeat {
    rep <- rep + 1
    prev.x <- cur.x
    delta.x <- abs(delta.x)
    order.x <- order(delta.x, decreasing=TRUE)
    dims <- partition2(nx, order.x, delta.x, subspc.min, subspc.max)

    ## Simplex loop.
    for ( i in seq_along(dims) ) {
      j <- dims[[i]]
      res <- tryCatch(simplex(constrain2(func, cur.x, j), cur.x[j],
                              step.size[j], max.eval-ne, y.init=cur.y,
                              dx=control$dx[j]),
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

    delta.x <- cur.x - prev.x
    done <- any(pmax(abs(delta.x), abs(step.size)*psi)/
                pmax(abs(cur.x),1) <= tolerance)
      
    if ( done ) { ## Tolerance is satisfied
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
  list(par=cur.x, value=cur.y, count=ne, convergence=flag, message=msg)
}

## Simplified simplex, following the code directly.
simplex <- function(func, x.init, scale, max.eval=1000, alpha=1.0,
                    beta=0.5, gamma=2.0, delta=0.5, psi=0.25,
                    y.init=NULL, dx=NULL) {
  simplexPrecisionError <- function(x, y) {
    e <- simpleError("Limit of machine precision reached")
    e$ret <- list(par=x, val=y, flag=1, ne=ne)
    class(e) <- c("simplexPrecisionError", class(e))
    stop(e)
  }
  newpt <- function(coef, x.base, x.old) {
    x.new <- x.base + coef * (x.base - x.old)
    ## Enforces grid.
    x.new[on.grid] <-
      round.to.grid(x.new[on.grid], x.init[on.grid], dx[on.grid])
    if ( identical(x.new, x.base) || identical(x.new, x.old) )
      simplexPrecisionError(xx[ind.lo,], yy[ind.lo])
    x.new
  }

  nx <- length(x.init)
  ne <- 0

  if ( is.null(dx) ) {
    on.grid <- rep(FALSE, nx)
  } else {
    on.grid <- !is.na(dx) & dx > 0
    scale[on.grid] <- pmax(round.to.grid(scale[on.grid], 0,
                                         dx[on.grid]), dx[on.grid])
  }

  if ( length(scale) != nx )
    stop("Invalid length scale")

  if ( is.null(y.init) ) {
    y.init <- func(x.init)
    ne <- 1
  }

  xx <- rbind(x.init, t(x.init + diag(nx) * scale), deparse.level=0)
  xl <- matrix.to.list(xx)
  if ( any(unlist(lapply(xl[-1], identical, x.init))) )
    simplexPrecisionError(x.init, y.init)
  yy <- c(y.init, unlist(lapply(xl[-1], func)))
  ne <- ne + nx

  ord <- order(yy)
  ind.lo  <- ord[1]
  ind.sec <- ord[nx]  
  ind.hi  <- ord[nx+1]

  tolerance <- psi^2 * sum((xx[ind.hi,]-xx[ind.lo,])^2)

  repeat {
    centroid <- colMeans(xx[-ind.hi,])

    ## "Reflection"
    x.ref <- newpt(alpha, centroid, xx[ind.hi,])
    y.ref <- func(x.ref)
    ne <- ne + 1

    if ( y.ref < yy[ind.lo] ) { # new best - keep pushing
      ## "Expansion"
      x.exp <- newpt(-gamma, centroid, x.ref)
      y.exp <- func(x.exp)
      ne <- ne + 1

      if ( y.exp < y.ref ) {
        xx[ind.hi,] <- x.exp
        yy[ind.hi] <- y.exp
      } else {
        xx[ind.hi,] <- x.ref
        yy[ind.hi] <- y.ref
      }
    } else if ( y.ref < yy[ind.sec] ) {
      ## Neither the best nor the worst (in the remaining simplex):
      ## accept current point
      xx[ind.hi,] <- x.ref
      yy[ind.hi] <- y.ref
    } else { ## Would be the current worst point.
      ## "Contraction"
      x.base <- if ( yy[ind.hi] < y.ref ) xx[ind.hi,] else x.ref
      x.con <- newpt(-beta, centroid, x.base)
      y.con <- func(x.con)
      ne <- ne + 1

      if ( y.con < min(y.ref, yy[ind.hi]) ) {
        ## Successful contraction!
        xx[ind.hi,] <- x.con
        yy[ind.hi] <- y.con
      } else {
        ## "Massive Contraction"
        for ( i in seq_len(nx + 1)[-ind.lo] ) {
          xx[i,] <- newpt(-delta, xx[ind.lo,], xx[i,])
          yy[i] <- func(xx[i,])
        }
        ne <- ne + nx
      }
    }

    ord <- order(yy)
    ind.lo <- ord[1]
    ind.sec <- ord[nx]
    ind.hi <- ord[nx+1]
    
    if ( ne > max.eval ) {
      flag <- -1
      break
    } else if ( sum((xx[ind.hi,]-xx[ind.lo,])^2) <= tolerance ) {
      flag <- 0
      break
    }
  }

  list(par=xx[ind.lo,], val=yy[ind.lo], flag=flag, ne=ne)
}

## This is a direct copy of the Matlab partition.m file's function
## partitionx().
## Description of what is going on:
##   "Partition the space defined by deltax[order] into a group of
##   sets, such that subs_max <= k <= subs_max for all sets of size k.
##   Where possible, pick combinations of points which maximize a-b,
##   where a is the norm of the new set /k and b is the norm of the
##   remainder over (nx-k)."
partition <- function(nx, order, dx, subs.min, subs.max) {
  num.subspaces <- 0
  nused <- 0
  nleft <- nx
  asleft <- sum(dx)
  subsp.dims <- rep(NA, nx)

  while ( nused < nx ) {
    num.subspaces <- num.subspaces + 1
    as1 <- sum(dx[order[(nused+1):(nused+subs.min-1)]])
    gapmax <- -1
    for ( ns1 in subs.min:(min(subs.max,nleft)) ) {
      as1 <- as1 + dx[order[nused+ns1]]
      ns2 <- nleft - ns1
      if ( ns2 > 0 ) {
        if ( ns2 >= (trunc((ns2-1)/subs.max+1)*subs.min) ) {
          as2 <- asleft - as1
          gap <- as1/ns1 - as2/ns2
          if ( gap > gapmax ) {
            gapmax <- gap
            subsp.dims[num.subspaces] <- ns1
            as1max <- as1
          }
        }
      } else {
        if ( as1/ns1 > gapmax ) {
          subsp.dims[num.subspaces] <- ns1
          return(list(num.subspaces, subsp.dims[1:num.subspaces]))
        }
      }
    }

    nused <- nused + subsp.dims[num.subspaces]
    nleft <- nx - nused
    asleft <- asleft - as1max
  }

  list(num.subspaces, subsp.dims[1:num.subspaces])
}

## This just wraps the above partition() function in a more R-ish way
## (at the cost of some speed due to the call to split()), which makes
## the book-keeping easier.
partition2 <- function(nx, order.x, delta.x, subspc.min, subspc.max) {
  tmp <- partition(nx, order.x, delta.x, subspc.min, subspc.max)
  ret <- split(order.x[seq_len(nx)], rep(seq_len(tmp[[1]]), tmp[[2]]))
  names(ret) <- NULL
  ret
}

## xmid may be any point that is on the grid.
round.to.grid <- function(x, xmid, dx)
  round((x - xmid)/dx)*dx + xmid
