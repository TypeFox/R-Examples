plot.fdSmooth <- function(x, y, Lfdobj=0, href=TRUE, titles=NULL,
                          xlim=NULL, ylim=NULL, xlab=NULL,
                          ylab=NULL, ask=FALSE, nx=NULL, axes=NULL,
                          ...){
  plot(x$fd, y, Lfdobj=Lfdobj, href=href, titles=titles,
       xlim=xlim, ylim=ylim, xlab=xlab,
       ylab=ylab, ask=ask, nx=nx, axes=axes, ...)
}

plot.fdPar <- function(x, y, Lfdobj=0, href=TRUE, titles=NULL,
                    xlim=NULL, ylim=NULL, xlab=NULL,
                    ylab=NULL, ask=FALSE, nx=NULL, axes=NULL, ...){
  plot(x$fd, y, Lfdobj=Lfdobj, href=href, titles=titles,
       xlim=xlim, ylim=ylim, xlab=xlab,
       ylab=ylab, ask=ask, nx=nx, axes=axes, ...)
}

plot.fd <- function(x, y, Lfdobj=0, href=TRUE, titles=NULL,
                    xlim=NULL, ylim=NULL, xlab=NULL,
                    ylab=NULL, ask=FALSE, nx=NULL, axes=NULL, ...)
{
#  -----------------------------------------------------------------------
#       plot for fd class
#  -----------------------------------------------------------------------

  #  Plot a functional data object fdobj.
  #  Arguments:
  #  fdobj     ... a functional data object
  #  Lfdobj    ... linear differental operator to be applied to fdobj before
  #             plotting
  #  HREF   ... If TRUE, a horizontal dotted line through 0 is plotted.
  #  The argument ASK, if TRUE, causes the curves to be displayed one at a time.
  #  NX     ... The number of sampling points to use for
  #             plotting.  (default 101)

  #  The remaining optional arguments are the same as those available
  #     in the regular "plot" function.

  #  Note that for multivariate fdobj, a suitable matrix of plots
  #    must be set up before calling plot by using something such as
  #    par(mfrow=c(1,nvar),pty="s")

# last modified July 6, 2012 by Spencer Graves
#   for argvals of class Date and POSIXct
# Prevously modified 2 May 2012 by Jim Ramsay
##
## 1.  Basic checks
##
  fdobj <- x
  if (!(inherits(fdobj, "fd"))) stop(
		"First argument is not a functional data object.")
  {
    if(is.null(axes)){
      if(is.null(fdobj$basis$axes)){
        Axes <- TRUE
        axFun <- FALSE
      }
      else{
        if(!inherits(fdobj$basis$axes, 'list'))
          stop('fdobj$basis$axes must be a list;  ',
               'class(fdobj$basis$axes) = ', class(fdobj$basis$axes))
        if(!(inherits(fdobj$basis$axes[[1]], 'character') ||
             inherits(fdobj$basis$axes[[1]], 'function') ) )
          stop('fdobj$basis$axes[[1]] must be either a function or the ',
               'name of a function;  class(fdobj$basis$axes[[1]]) = ',
               class(fdobj$basis$axes[[1]]) )
        Axes <- FALSE
        axFun <- TRUE
        axList <- c(fdobj$basis$axes, ...)
      }
    }
    else{
      if(is.logical(axes)){
        Axes <- axes
        axFun <- FALSE
      }
      else{
        if(!inherits(axes, 'list'))
          stop('axes must be a logical or a list;  class(axes) = ',
               class(axes))
        if(!(inherits(axes[[1]], 'character') ||
             inherits(axes[[1]], 'function') ) )
          stop('axes[[1]] must be either a function or the ',
               'name of a function;  class(axes[[1]]) = ',
               class(axes[[1]]) )
        Axes <- FALSE
        axFun <- TRUE
        axList <- c(axes, ...)
      }
    }
  }
  #  check fdobj

  #  check Lfdobj

  Lfdobj <- int2Lfd(Lfdobj)
  if (!inherits(Lfdobj, "Lfd")) stop(
      "Second argument is not a linear differential operator.")

  #  extract dimension information

  coef   <- fdobj$coefs
  coefd  <- dim(coef)
  ndim   <- length(coefd)
  # Number of basis functions
  nbasis    <- coefd[1]
  if(is.null(nx)) nx <- max(c(501,10*nbasis + 1))
  # Number of functional observations
  nrep   <- coefd[2]
  if (ndim > 2) nvar <- coefd[3] else nvar <- 1

  #  get basis information

  basisobj <- fdobj$basis
  rangex   <- basisobj$rangeval
  #  set up a set of argument values for the plot

  if (missing(y)) {
    y <- nx
  } else {
    if(is.numeric(y)) y <- as.vector(y)
  }

  Y <- y
  if (length(y) == 1) {
    if (y >= 1) {
      y <- seq(rangex[1],rangex[2],len=round(y))
    } else {
      stop("'y' a single number less than one.")
    }
  }
  if (min(y) < rangex[1] || max(y) > rangex[2])
    stop("Values in Y are outside the basis range.")
  if (is.null(xlim)){
      xlim <- rangex
  } else {
      rangex[1] <- max(rangex[1], xlim[1])
      rangex[2] <- min(rangex[2], xlim[2])
      if(length(Y)==1)
          y <- seq(rangex[1],rangex[2],len=round(Y))
  }

  #  evaluate LFDOBJ(FDOBJ) at the argument values

  fdmat    <- eval.fd(y, fdobj, Lfdobj)
  rangey   <- range(fdmat)
  if (is.null(ylim)) ylim <- rangey

  #  set up axis labels and,
  #  optionally, caselabels and variable labels

  fdnames      = fdobj$fdnames
  fdlabelslist = fdlabels(fdnames, nrep, nvar)

# Ramsay 2008.08.26
  xlabel    = fdlabelslist$xlabel
  ylabel    = fdlabelslist$ylabel
  casenames = fdlabelslist$casenames
  varnames  = fdlabelslist$varnames

# Graves 2008.07.04
#  xlabel   <- names(fdobj$fdnames)[[1]]
#  ylabel   <- names(fdobj$fdnames)[[3]]
#  if (is.character(xlabel) == FALSE) xlabel <- ""
#  if (is.character(ylabel) == FALSE) ylabel <- ""

  #  check xlab and ylab
  if (is.null(xlab)) xlab <- xlabel
  if (is.null(ylab)) ylab <- ylabel
#  if (missing(xlab)) xlab <- xlabel
#  if (missing(ylab)) ylab <- ylabel
#  crvnames <- fdobj$fdnames[[2]]
#  varnames <- fdobj$fdnames[[3]]
# Don't ask for the first plot, but do for later plots if(ask)
#  op <- par(ask=FALSE)
# Don't ask for the first plot,
# but if ask==TRUE, do ask for succeeding plots
#  on.exit(par(op))
# A single line?

  # A single line?
  if (ndim < 2) {
    plot (y, fdmat, type="l", xlim=xlim, ylim=ylim,
          xlab=xlab, ylab=ylab, axes=Axes, ...)
    if(axFun)
        do.call(axList[[1]], axList[-1])
#   Ramsay 2008.08.26
    if (zerofind(fdmat) && href) abline(h=0,lty=2)
#   Graves 2008.07.04
#    if (zerofind(ylim) && href) abline(h=0,lty=2)
  }
  # Several copies of one function?
  if (ndim ==2 ) {
    if (!ask) {
      matplot(y, fdmat, type="l",
              xlim=xlim,   ylim=ylim,
              xlab=xlab, ylab=ylab, axes=Axes, ...)
      if(axFun)
          do.call(axList[[1]], axList[-1])
#   Ramsay 2008.08.26
      if (zerofind(fdmat) && href) abline(h=0,lty=2)
#   Graves 2008.07.04
#    if (zerofind(ylim) && href) abline(h=0,lty=2)
    } else  {
#   Graves 2008.07.04:  par, cat absent from Ramsay 2008.08.26
      op <- par(ask=FALSE)
# Don't ask for the first plot,
# but if ask==TRUE, do ask for succeeding plots
      on.exit(par(op))
      cat('Multiple plots:  Click in the plot to advance to the next')
#      op <- par(ask = TRUE)
#      on.exit(par(op))
      for (irep in 1:nrep) {
        plot (y, fdmat[,irep], type="l",
              xlim=xlim, ylim=ylim,
              xlab=xlab, ylab=ylab, axes=Axes, ...)
        if(axFun)
            do.call(axList[[1]], axList[-1])
        if(irep<2) par(ask=ask)

#        if (zerofind(ylim) && href) abline(h=0,lty=2)
#        if (!is.null(titles)) title(titles[irep])
#        else title(paste(crvnames[irep]))

        if (!is.null(casenames)) title(casenames[irep])
        else                     title(paste("Case",irep))
#        if (zerofind(fdmat[,irep]) && href) abline(h=0,lty=2)
        if (zerofind(ylim) && href) abline(h=0,lty=2)

#        mtext("Click in graph to see next plot", side=3, outer=FALSE)
#        text("",locator(1))
      }
    }
  }
  # Possibly multiple copies of different functions
  if (ndim == 3) {
    if (!ask) {
      for (ivar in 1:nvar) {
        matplot (y, fdmat[,,ivar], type="l",
                 xlim=xlim, ylim=ylim,
                 xlab=xlab, ylab=ylab, ask=FALSE, axes=Axes, ...)
        if(axFun)
            do.call(axList[[1]], axList[-1])
        if (!is.null(varnames)) title(varnames[ivar])
        else                    title(paste("Variable",ivar))
#        if (zerofind(fdmat[,,ivar]) && href) abline(h=0,lty=2)
        if (zerofind(ylim) && href) abline(h=0,lty=2)
      }
    } else {
      op <- par(ask=FALSE)
# Don't ask for the first plot,
# but if ask==TRUE, do ask for succeeding plots
      on.exit(par(op))
      cat('Multiple plots:  Click in the plot to advance to the next')

      for (irep in 1:nrep) {
        for (ivar in 1:nvar) {
          plot(y,fdmat[,irep,ivar],type="l",
               xlim=xlim, ylim=ylim,
               xlab=xlab, ylab=ylab, axes=Axes, ...)
          if(axFun)
              do.call(axList[[1]], axList[-1])
          if (!is.null(casenames)) titlestr = casenames[irep]
          else                     titlestr = paste("Case",irep)
          if (!is.null(varnames)) {
             titlestr = paste(titlestr,"  ",varnames[ivar])
          } else {
             titlestr = paste(titlestr,"  ","Variable",ivar)
          }
          title(titlestr)
#          if (zerofind(fdmat[,irep,ivar]) && href) abline(h=0,lty=2)
          if (zerofind(ylim) && href) abline(h=0,lty=2)
#          if (!is.null(titles)) title(titles[irep])
#          else title(paste("Curve", irep, varnames[ivar]))

#          mtext("Click in graph to see next plot", side=3, outer=FALSE)
#          text("",locator(1))
        }
      }
    }
  }
#  invisible(NULL)
# This used to return 'invisible(NULL)'.
# However, with R 2.7.0 under XEmacs with ESS,
# it sometimes failed to plot.  I changed the return value,
# and the problem disappeared.
  'done'
}

#  --------------------------------------------------------------------

zerofind <- function(fmat)
{
  frng <- range(fmat)
  if (frng[1] <= 0 && frng[2] >= 0) zeroin <- TRUE else zeroin <- FALSE
  return(zeroin)
}



