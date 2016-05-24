#' @title Various utility functions.
#'
#' @description \emph{The various utility functions are:}
#'
#' @author A.J. Barbour
#' @name psd-utilities
#' @seealso \code{\link{psd-package}}, \code{\link{as.tapers}}, \code{\link{modulo_floor}}
#' 
#' @example inst/Examples/rdex_utilities.R
NULL

#' @description \code{\link{vardiff}} returns the variance of the first (or second) 
#' difference of the series. \code{\link{varddiff}} is a convenience wrapper
#' to return variance for the second difference.
#' @rdname psd-utilities
#' @export
#' @param double.diff logical; should the double difference be used instead?
vardiff <- function(x, double.diff=FALSE){
  dorder <- ifelse(double.diff, 2, 1)
  stats::var(diff(x, differences=dorder))
}
#' @rdname psd-utilities
#' @export
varddiff <- function(x) UseMethod('varddiff')
#' @rdname psd-utilities
#' @export
varddiff.spec <- function(x) varddiff(as.vector(x[['spec']]))
#' @rdname psd-utilities
#' @export
varddiff.default <- function(x) vardiff(x, double.diff=TRUE)

#' @rdname psd-utilities
#' @description \code{\link{create_poly}} generates an x-y sequence compatible for use with \code{\link{polygon}}
#' @param x,y objects; in \code{\link{create_poly}} these are the vectors used to
#' create a \code{\link{polygon}}-compatible sequence (\code{x} is sorted by default);
#' in \code{\link{mod}} these are the "numerator" and "denominator", respectively.
#' @param dy numeric; the distance from \code{y} to the top and bottom of
#' the polygonal surfaces; see \code{from.lower}
#' @param from.lower logical; should the bottom be \code{y} instead of \code{y+dy}, so that
#' \code{dy} represents the distance from the lower surface?
#' @export
create_poly <- function(x, y, dy, from.lower=FALSE){
  x <- sort(x)
  xx <- c(x, rev(x))
  yy <- if (from.lower){
    c(y, rev(y+dy))
  } else {
    c(y+dy, rev(y-dy))
  }
  return(data.frame(x.x=xx, y.y=yy))
}

#' @description \code{\link{dB}} returns an object converted to decibels.
#' @details 
#' Decibels are defined as \eqn{10 \log{}_{10} \frac{X_1}{X_2}}, 
#' unless \code{is.power=TRUE} in which \eqn{\mathrm{db} X^2 \equiv 20 \log{}_{10} X^2}
#' @rdname psd-utilities
#' @param Rat numeric; the values -- ratios -- to convert to decibels (\code{dB}).
#' @param invert logical; assumes \code{Rat} is already in decibels, so return ratio
#' @param pos.only logical; if \code{invert=FALSE}, sets negative or zero values to NA
#' @param is.power logical; should the factor of 2 be included in the decibel calculation?
#' @export
#' @aliases decibels db
dB <- function(Rat, invert=FALSE, pos.only=TRUE, is.power=FALSE){
  CC <- ifelse(is.power, 20, 10)
  if (invert) {
    10 ** (Rat/CC)
  } else {
    if (pos.only) Rat[Rat <= 0] <- NA
    CC * log10(Rat)
  }
}

#' @description \code{\link{vector_reshape}} reshapes a vector into another vector.
#' @rdname psd-utilities
#' @name vector_reshape
#' @param vec.shape  choice between horizontally-long or vertically-long vector.
#' @return \code{vector_reshape} returns a "reshaped" vector, meaning it has
#' had it's dimensions changes so that it has either one row 
#' (if \code{vec.shape=="horizontal"}), or one column (\code{"vertical"}).
#' @export
vector_reshape <- function(x, vec.shape=c("horizontal","vertical")) {
  x <- as.vector(x)
  vec.shape <- match.arg(vec.shape)
  nrow <- switch(vec.shape, "horizontal"=1, "vertical"=length(x))
  return(matrix(x, nrow=nrow))
}

#' @description \code{\link{colvec}} returns the object as a vertically long vector; whereas
#' \code{\link{rowvec}} returns the object as a horizontally long vector.
#' @details \code{colvec, rowvec} are simple wrapper functions to \code{vector_reshape}.
#' @rdname psd-utilities
#' @export
colvec <- function(x) vector_reshape(x, "vertical")

#' @rdname psd-utilities
#' @export
rowvec <- function(x) vector_reshape(x, "horizontal")

#' @description \code{\link{is.spec}} and \code{\link{is.amt}} report whether an object has class \code{'spec'} or \code{'amt'}, as
#' would one returned by, for example, \code{\link{spectrum}} or \code{\link{psdcore}}.
#' 
#' \code{\link{is.tapers}} reports whether an object has class \code{'tapers'}, as
#' would one returned by, for example, \code{\link{as.tapers}}.
#' 
#' @param Obj  An object to test for class inheritance.
#' @return \code{is.spec}, \code{is.amt}, and \code{is.tapers} return the output of \code{\link{inherits}}.
#' @rdname psd-utilities
#' @export
is.spec <- function(Obj) inherits(Obj, "spec")

#' @rdname psd-utilities
#' @export
is.amt <- function(Obj) inherits(Obj, 'amt')

#' @rdname psd-utilities
#' @export
is.tapers <- function(Obj) inherits(Obj, "tapers")

#' @description \code{\link{na_mat}} populates a matrix of specified dimensions 
#' with \code{NA} values.
#' @rdname psd-utilities
#' @param nrow,ncol integer; the number of rows and/or columns to create
#' @return \code{na_mat} returns a matrix of dimensions \code{(nrow,ncol)} with
#' \code{NA} values, the representation of which is set by \code{NA_real_}
#' @export
na_mat <- function(nrow, ncol=1) {matrix(NA_real_, nrow, ncol)}

#' @description \code{\link{zeros}} populate a column-wise matrix with zeros; whereas,
#' \code{\link{ones}} populates a column-wise matrix with ones.  \emph{Note that 
#' \code{n} is enforced to be at least 1 for both functions.}
#' @rdname psd-utilities
#' @export 
zeros <- function(nrow) {
  nrow <- max(1., abs(nrow))
  matrix(0., nrow=nrow, ncol=1)
}

#' @rdname psd-utilities
#' @export
ones <- function(nrow) {
  nrow <- max(1., abs(nrow))
  matrix(1., nrow=nrow, ncol=1)
}

#' @description \code{\link{mod}} finds the modulo division of two values
#' 
#' @details Modulo division has higher order-of-operations ranking than other
#' arithmetic operations; hence, \code{x + 1 \%\% y} is equivalent to
#' \code{x + (1 \%\% y)} which can produce confusing results. \code{mod}
#' is simply a series of \code{trunc} commands which
#' reduces the chance for unintentionally erroneous results.
#' 
#' @note The performance of \code{\link{mod}} has not been tested against the 
#' \code{\%\%} arithmetic method -- it may or may not be slower for large
#' numeric vectors.
#' 
#' @references For \code{\link{mod}}: see Peter Dalgaard's explanation of 
#' the non-bug (#14771) I raised (instead I should've asked it on R-help): 
#' \url{https://bugs.r-project.org/bugzilla3/show_bug.cgi?id=14771\#c2}
#' 
#' @rdname psd-utilities
#' @export
#' @return \code{mod} returns the result of a modulo division, which is 
#' equivalent to \code{(x) \%\% (y)}.
mod <- function(x, y) {
  stopifnot(is.numeric(c(x, y)))
  ## modulo division
  x1 <- trunc( trunc(x/y) * y)
  z <- trunc(x) - x1
  return(z)
}


#' Numerical derivatives of a series based on its smooth-spline representation
#' 
#' @description This computes the numerical derivatives of a spline 
#' representation of the input series; differentiation of spline curves is 
#' numerically efficient.
#' 
#' @details
#' With smoothing, the numerical instability for "noisy" data can be drastically
#' reduced, since spline curves are inherently (at least) twice differentiable. 
#' 
#' @author A.J. Barbour
#' @name splineGrad
#' @param dseq  numeric; a vector of positions for \code{dsig}.
#' @param dsig  numeric; a vector of values (which will have a spline fit to them).
#' @param plot.derivs  logical; should the derivatives be plotted?
#' @param ... additional arguments passed to \code{\link{smooth.spline}}
#' @return A matrix with columns representing \eqn{x, f(x), f'(x), f''(x)}
#' @export
#' @seealso \code{\link{smooth.spline}}, \code{\link{constrain_tapers}}
#' @example inst/Examples/rdex_splinegrad.R
splineGrad <- function(dseq, dsig, ...) UseMethod("splineGrad")

#' @rdname splineGrad
#' @aliases splineGrad.default
#' @export
splineGrad.default <- function(dseq, dsig, plot.derivs=FALSE, ...){
  #
  # Use spline interpolation to help find an emprirical gradient
  # (reduces numerical instability)
  #
  # @dseq: the sequence (index) for @dsig, the signal
  # output is the same length as the input
  #
  # create a weighted cubic spline
  smspl <- stats::smooth.spline(dseq, dsig, ...)
  # and a function
  SPLFUN <- stats::splinefun(smspl$x, smspl$y)
  # ?splinefun:
  # splinefun returns a function with formal arguments x and deriv, 
  # the latter defaulting to zero. This function can be used to 
  # evaluate the interpolating cubic spline (deriv = 0), or its 
  # derivatives (deriv = 1, 2, 3) at the points x, where the spline 
  # function interpolates the data points originally specified. This 
  # is often more useful than spline.
  #
  #   seq.rng <- range(dseq)
  #   from <- seq.rng[1]
  #   to <- seq.rng[2]
  #   n <- length(dseq)
  #
  # signal spline
  #fsig <<- SPLFUN(dseq)
  #   FD0 <- function(){graphics::curve(SPLFUN(x), from=from, to=to, n=n, add=NA)}
  #   fsig <<- FD0()
  #
  # first deriv
  #   FD1 <- function(){graphics::curve(SPLFUN(x,deriv=1), from=from, to=to, n=n, add=NA)}
  #   fsigderiv <<- FD1()
  fsigderiv <- SPLFUN(dseq, deriv=1)
  #
  # second deriv
  #   FD2 <- function(){graphics::curve(SPLFUN(x,deriv=2), from=from, to=to, n=n, add=NA)}
  #   fsigderiv2 <<- FD2()
  fsigderiv2 <- SPLFUN(dseq, deriv=2)
  # how does the first deriv of the first deriv compare to the second?
  #   smspl.alt <- stats::smooth.spline(dseq, fsigderiv, ...)
  #   SPLFUN.alt <- stats::splinefun(smspl.alt$x, smspl.alt$y)
  #   fsigderiv2.alt <<- SPLFUN.alt(dseq, deriv=1)
  #   print(all.equal(fsigderiv2,fsigderiv2.alt))
  # [1] TRUE
  #   plot(fsigderiv2,fsigderiv2.alt, asp=1)
  ##
  toret <- data.frame(x=dseq, y=dsig, 
                      dydx=fsigderiv, 
                      d2yd2x=fsigderiv2)
  #d2yd2x.alt=fsigderiv2.alt)
  #
  if (plot.derivs){
    #     yl.u <- max(c(dsig,fsigderiv,fsigderiv2))#,fsigderiv2.alt))
    #     yl.l <- min(c(dsig,fsigderiv,fsigderiv2))#,fsigderiv2.alt))
    nr <- 3 # f, f', f''
    mar.multi <- c(2., 5.1, 2, 2.1)
    oma.multi <- c(6, 0, 5, 0)
    oldpar <- par(mar = mar.multi, oma = oma.multi, mfcol = c(nr, 1))
    on.exit(par(oldpar))
    par(las=1)
    plot(dseq, dsig, cex=0.6, pch=3,
         #ylim=1.1*c(yl.l,yl.u),
         xaxs="i", yaxs="i",
         xlab="x", ylab="f(x)",
         main=sprintf("splineGrad: signal and weighted cubic-spline fit (spar = %g)",smspl$spar))
    lines(y ~ x, toret, col="dark grey", lwd=1.5)
    plot(dydx ~ x, toret, 
         xaxs="i", yaxs="i",
         main="first derivative",
         col="red", type="s", lwd=2.4)
    plot(d2yd2x ~ x, toret, 
         xaxs="i", yaxs="i",
         main="second derivative", xlab="x",
         col="blue", type="s", lwd=2.4, lty=3)
  }
  return(invisible(toret))
}
