plot.basisfd <- function(x, knots=TRUE, axes=NULL, ...) {
  basisobj <- x
#  plot a basis object

# last modified July 3, 2012 by Spencer Graves for Date and POSIXct
# previously modified 26 May 2012 by Jim Ramsay
##
##  check BASISOBJ
##
  if (!inherits(basisobj, "basisfd"))
    stop("argument x is not a basis object.")
##
## Check dots
##
  dot.args <- list(...)
  {
    if(is.null(axes)){
      if(is.null(x$axes)){
        dot.args$axes <- TRUE
        axFun <- FALSE
      } else {
        if(!inherits(x$axes, 'list'))
          stop('x$axes must be a list;  class(x$axes) = ',
               class(x$axes))
        if(!(inherits(x$axes[[1]], 'character') ||
             inherits(x$axes[[1]], 'function') ) )
          stop('x$axes[[1]] must be either a function or the ',
               'name of a function;  class(x$axes[[1]]) = ',
               class(x$axes[[1]]) )
        axList <- c(x$axes, ...)
        dot.args$axes <- FALSE
        axFun <- TRUE
      }
    } else{
      if(is.logical(axes)){
        dot.args$axes <- axes
        axFun <- FALSE
      } else {
        if(!inherits(axes, 'list'))
          stop('axes must be a logical or a list;  class(axes) = ',
               class(axes))
        if(!(inherits(axes[[1]], 'character') ||
             inherits(axes[[1]], 'function') ) )
          stop('axes[[1]] must be either a function or the ',
               'name of a function;  class(axes[[1]]) = ',
               class(axes[[1]]) )
        axList <- c(axes, ...)
        dot.args$axes <- FALSE
        axFun <- TRUE
      }
    }
  }

  if(is.null(dot.args$xlab))dot.args$xlab <- ''
  if(is.null(dot.args$ylab))dot.args$ylab <- ''
##
## get basis info
##
  nbasis   <- basisobj$nbasis

  if(is.null(dot.args$type))dot.args$type <- 'l'

  if(is.null(dot.args$lty))
    dot.args$lty <- rep(1:3, max(1, nbasis/3))

  nx       <- max(501,10*nbasis)

  {
    if(is.null(dot.args$xlim)){
        rangex <- basisobj$rangeval
    } else {
      rangex <- dot.args$xlim
# shrink rangex so it does not exceed basisobj$rangeval
      rangex[1] <- max(basisobj$rangeval[1], dot.args$xlim[1])
      rangex[2] <- min(basisobj$rangeval[2], dot.args$xlim[2])
    }
  }

  argvals  <- seq(rangex[1],rangex[2],len=nx)
  basismat <- eval.basis(argvals, basisobj)

#  minval   <- min(basismat)
#  maxval   <- max(basismat)
#  if (minval == maxval) {
#    if (abs(minval) < 1e-1) {
#      minval <- minval - 0.05
#      maxval <- maxval + 0.05
#    } else {
#      minval <- minval - 0.05*minval
#      maxval <- maxval + 0.05*minval
#    }
#  }

#  matplot (argvals, basismat, type="l", lty=ltype,
#           xlab=xlabel, ylab=ylabel, cex=cexval,
#           xlim=c(argvals[1],argvals[nx]),
#           ylim=, ...)

  dot.args$x <- range(argvals)
  dot.args$y <- range(basismat)
  type <- dot.args$type
  dot.args$type <- 'n'

#  do.call('matplot', dot.args)
  do.call(plot, dot.args)
#
  dot.args$x <- argvals
  dot.args$y <- basismat
  dot.args$type <- type
  do.call(matlines, dot.args)
# knots?
  if(knots && (x$type=='bspline'))
    abline(v=knots(x), lty='dotted', col='red')
# axes?
  if(axFun)
    do.call(axList[[1]], axList[-1])
}
