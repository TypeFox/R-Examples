register.fd0 <- function(y0fd, yfd=NULL, ...)
{
##
## 1.  Check y0fd and yfd
##
#  check classes of first two arguments
  if(missing(y0fd)) stop('y0fd missing with no default.')
  if(is.null(yfd)){
    yfd  <- y0fd
    y0fd <- NULL
  }
#
  if (!(inherits(yfd, "fd")))
      stop("'yfd' must be a functional data object.  ",
           "Instead, class(yfd) = ", class(yfd))
#
  if(is.null(y0fd)) {
    y0fd <- mean.fd(yfd)
  } else {
    if (!(inherits(y0fd, "fd")))
        stop("First argument is not a functional data object.",
             "Instead, class(y0fd) = ", class(y0fd))
  }

#  check functions to be registered

  ydim   <- dim(yfd$coefs)
  ncurve <- ydim[2]
  ndimy  <- length(ydim)
  if (ndimy == 3) {
      nvar <- ydim[3]
  } else {
      nvar <- 1
  }

  if (ndimy > 3) stop("'yfd' is more than 3-dimensional.")

#  Check target function

  y0coefs0 <- y0fd$coefs
  y0dim0   <- dim(y0coefs0)
  ndimy00  <- length(y0dim0)
  if (ndimy00 > ndimy) stop("Y0FD has more dimensions than YFD")
#  Determine whether the target function is full or not
#  if (y0dim0[2] == 1) {
#      fulltarg <- FALSE
#  } else {
#      if (y0dim0[2] == ydim[2]) {
#          fulltarg <- TRUE
#      } else {
#          stop("Second argument of Y0FD not correct.")
#      }
#  }
  if (ndimy00 == 3 && ydim[3] != y0dim0[3]) stop(
      "Third dimension of YOFD does not match that of YFD.")

  y0rangeval <- y0fd$basis$rangeval
  if(any(yfd$basis$rangeval != y0rangeval))
    stop('y0fd and yfd do not have the same rangeval')
##
## 2.  Internal functions for optimize
##
#  2.1.  integrand function
  dy2 <- function(x, x0, y, y0){
#  2.1.1.  domain of y, y0 adjusted by x0
    y0range <- y0$basis$rangeval
    xlim <- y0range
    if(x0>0){
      xlim[1] <- xlim[1]+x0
    } else xlim[2] <- xlim[2]+x0
#
    d.x <- diff(xlim)
    if(d.x<=0){
#      warning('range reduced to 0 in ss.dy2;  returning Inf')
      return(Inf)
    }
#  2.1.2.  rescale x
    x. <- xlim[1]+d.x*x
    y. <- eval.fd(x.-x0, y)
    y0. <- eval.fd(x., y0)
    (y.-y0.)^2 / d.x
  }
#  2.2.  objective function
  ss.dy2 <- function(x0, y, y0, ...){
    int <- integrate(dy2, 0, 1, x0=x0, y=y, y0=y0, ...)
    int$value
  }
##
## 3.  optimize
##
  x0 <- rep(NA, ncurve)
  dxlim <- diff(y0rangeval)
  for(ic in 1:ncurve){
    fiti <- optimize(ss.dy2, c(-dxlim, dxlim), y=yfd[ic], y0=y0fd, ...)
    x0[ic] <- fiti$minimum
  }
##
## 4.  Compute regfd = the registerd yfd, starting with basis
##
  b <- yfd$basis
#  4.1.  rangeval
  newrangeval <- c(b$rangeval[1]+max(0, x0),
                   b$rangeval[2]+min(0, x0) )
  if(diff(newrangeval)<=0){
    warning('the adjusted rangeval is reduced to nothing:  ',
            paste(newrangeval, collapse=', ') )
    regfd <- rep(NA, ncurve)
    for(ic in 1:ncurve){
      regfd[ic] <- eval.fd(newrangeval[1]-x0[ic], yfd[ic])
    }
    dregfd <- regfd-eval.fd(newrangeval[1], y0fd)
  } else {
    npts <- 300
#  4.2.  if B-spline:  new basis
    if(b$type=='bspline'){
      knotsgood <- ((newrangeval[1] < b$params) &
                    (b$params < newrangeval[2]) )
      goodknots <- c(newrangeval[1], b$params[knotsgood],
                     newrangeval[2])
      b <- create.bspline.basis(newrangeval, norder=norder(b),
                                   breaks=goodknots)
      npts <- npts+length(goodknots)
    }
#  4.3.  Fit each curve to this new basis
    xx <- seq(newrangeval[1], newrangeval[2], length=npts)
    yxx <- matrix(NA, npts, ncurve)
    for(ic in 1:ncurve)
      yxx[, ic] <- eval.fd(xx-x0[ic], yfd[ic])
#
    regfd <- Data2fd(xx, yxx, b)
##
## 5.  dregfd
##
    y0xx <- as.numeric(eval.fd(xx, y0fd))
    dyxx <- (yxx-y0xx)
    dregfd <- Data2fd(xx, dyxx, b)
  }
##
## 6.  Done
##
  rfd <- list(regfd=regfd, dregfd=dregfd, offset=x0)
  class(rfd) <- 'register.fd0'
  rfd
}

