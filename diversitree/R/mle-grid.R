## Grid-limited simplex.  This is not the same as the the Burmen et al
## "Grid restrained Nelder-Mead Algorithm" (doi:
## 10.1007/s10589-005-3912-z), but hopefully will be useful for
## problems where there is a natural "scale" that parameters should
## not be evaluated within.

simplex.grid <- function(func, x.init, scale, dx, max.eval=1000,
                         alpha=1.0, beta=0.5, gamma=2.0, delta=0.5,
                         psi=0.25, y.init=NULL) {
  simplexPrecisionError <- function(x, y) {
    e <- simpleError("Limit of machine precision reached")
    e$ret <- list(par=x, val=y, flag=1, ne=ne)
    class(e) <- c("simplexPrecisionError", class(e))
    stop(e)
  }
  newpt <- function(coef, x.base, x.old) {
    x.new <- x.base + coef * (x.base - x.old)
    ## Enforces grid.
    x.new[on.grid] <- x.old[on.grid] +
      round((x.new - x.init)[on.grid]/dx[on.grid])
    if ( identical(x.new, x.base) || identical(x.new, x.old) )
      simplexPrecisionError(xx[ind.lo,], yy[ind.lo])
    x.new
  }

  on.grid <- !is.na(dx) & dx > 0

  nx <- length(x.init)
  ne <- 0

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
                         
                         
