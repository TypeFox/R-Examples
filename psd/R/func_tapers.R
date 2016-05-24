#' Coerce an object into a \code{'tapers'} object.
#'
#' In a tapered spectrum estimation algorithm, it is
#' necessary to enforce rules on the number of tapers
#' that may be applied.
#'
#' Formal requirements enforced by this function are:
#' \itemize{
#' \item Non-zero.
#' \item Integer values.
#' \item Fewer than the half-length of the spectrum.
#' }
#' For example, we cannot apply
#' zero tapers (the result would be a raw periodogram)
#' or one million tapers (that would be absurd, and
#' violate orthogonality
#' conditions for any series less than two million terms long!).
#' 
#' An object with S3 class \code{'tapers'} is created;
#' this will have
#' a minimum number of tapers in each position
#' set by \code{min_taper}, and
#' a maximum number of tapers in each position
#' set by \code{max_taper}.
#' If \code{minspan=TRUE}, the bounded taper is fed through \code{\link{minspan}}
#' which will restrict the maximum tapers to less than or equal to
#' the half-length of the spectrum.
#'
#' Various classes can be coerced into a \code{'tapers'} object; those
#' tested sofar include: scalar, vector, matrix, data.frame, 
#' and list.  
#'
#' Multiple objects are concatenated into a single
#' vector dimension.  
#' 
#' Enabling \code{setspan} will only override
#' \code{max_taper} should it be larger than the half-width of the series.
#'
# @section Example of For example, if the object is 
# \code{list(x=c(1,2),y=c(3,4,5,0,1.1))} then the corresponding \code{'tapers'}
# objects for the following arguments are:
#
# \describe{
# \item{\emph{defaults}}{\code{[1,2,3,4,5,1,1]}}
# \item{\code{setspan=TRUE}}{\code{[1,2,3,3,3,1,1]}}
# \item{\code{max_taper=5}}{\code{[1,2,3,4,5,1,1]}}
# \item{\code{max_taper=5,setspan=TRUE}}{\code{[1,2,3,3,3,1,1]}}
# }
#'
#' @note No support (yet) for use of \code{min_taper,max_taper} as vectors, although
#' this could be quite desirable.
#'
#' @param x An object to set
#' @param min_taper Set all values less than this to this.
#' @param max_taper Set all values greater than this to this.
#' @param setspan logical; should the tapers object be passed through \code{\link{minspan}} before being returned?
#' @param record.last logical; should the \code{x} be saved to the \code{\link{psd-environment}} before coercion?
#' @export
#' @return An object with class \code{'taper'}
#' @author A.J. Barbour
#' @seealso \code{\link{is.tapers}}
#' @example inst/Examples/rdex_tapers.R
as.tapers <- function(x, min_taper=1, max_taper=NULL, setspan=FALSE, record.last=FALSE){
  # taper should be non-zero integer, since it represents the
  # number of tapered sections to average; hence, floor.
  x <- as.vector(unlist(x))
  stopifnot(!is.character(x))
  record <- 'last_as_tapers'
  if (record.last) psd_envAssign(record, x)
  
  # TODO: set na.rm until we are sure resample_fft_rcpp returns only finite values
  if (is.null(max_taper)) max_taper <- floor(max(x, na.rm=TRUE))
  if (!(min_taper*max_taper >= 1  &  max_taper >= min_taper)){
    stop('Bad taper limits:\tmin ', min_taper, "\tmax ", max_taper)
  }
  #
  x <- as.integer( pmin.int(max_taper, pmax.int(min_taper, round(x), na.rm=TRUE), na.rm=TRUE) )
  attr(x, 'last_recorded') <- ifelse(record.last, record, NA)
  attr(x, "n_taper_limits_orig") <- c(min_taper, max_taper)
  attr(x, "taper_positions") <- NA
  #
  if (setspan) x <- as.integer( minspan(x) )
  attr(x, "span_was_set") <- setspan
  attr(x, "n_taper_limits") <- range(x, na.rm=TRUE)
  #
  class(x) <- "tapers"
  #
  return(x)
}
#' @rdname as.tapers
#' @name tapers
#' @export
tapers <- as.tapers

###
###  Generic methods
###

#' @title Generic methods for objects with class \code{'tapers'}
#' @name tapers-methods
#' @author A.J. Barbour
#' @rdname tapers-methods
#' @docType methods
#'
#' @seealso \code{\link{as.tapers}}, \code{\link{constrain_tapers}}
#' @param x tapers object
#' @param xi optional vector for indices of \code{x}
#' @param object tapers object
#' @param lwd line width (default is 1.8)
#' @param col color of line (default is "red")
#' @param pch point character (default is "_")
#' @param cex point size (default is 1)
#' @param color.pal color palette to use (choices are: "Blues","Spectral")
#' @param ylim optional limits for y-axis
#' @param hv.lines logical; should horizontal (log2) and vertical reference lines be plotted?
#' @param log.y logical; should the vertical scale be logarithmic?
#' @param ... optional arguments
#' @return \code{plot} returns a list with names: \code{line.colors} (hex values)
#' @examples
#' ##
#' tap <- as.tapers(c(1:49,50:0)+rnorm(1e2))
#' print(tap)
#' print(summary(tap))
#' plot(tap)
#' # no arithmetic methods
#' tap <- as.tapers(tap/2)
#' lines(tap)
NULL

#' @rdname tapers-methods
#' @export
as.data.frame.tapers <- function(x, ...){
  df <- as.data.frame.numeric(x)
  names(df) <- "n.tapers"
  return(df)
}
#' @rdname tapers-methods
#' @export
data.frame.tapers <- as.data.frame.tapers

#' @rdname tapers-methods
#' @aliases print.tapers
#' @export
print.tapers <- function(x, ...){
  stopifnot(is.tapers(x))
  xh <- paste(as.character(head(x)), collapse=" ")
  xt <- paste(as.character(tail(x)), collapse=" ")
  cat(sprintf("'tapers' object: num. tapers applied by index\n\thead:  %s\n\t\t...\n\ttail:  %s\n",xh,xt))
  spans <- attr(x, "span_was_set")
  if (spans) message('span-set: TRUE')
}

#' @rdname tapers-methods
#' @aliases summary.tapers
#' @export
summary.tapers <- function(object, ...){
  stopifnot(is.tapers(object))
  toret <- summary.default(object)
  class(toret) <- "summary.tapers"
  return(toret)
}

#' @rdname tapers-methods
#' @aliases print.summary.tapers
#' @export
print.summary.tapers <- function(x, ...){
  cat("summary of tapers:\n")
  print(summary(x))
}

#' @rdname tapers-methods
#' @aliases lines.tapers
#' @export
lines.tapers <- function(x, lwd=1.8, col="red", ...){
  stopifnot(is.tapers(x))
  nt <- length(x)
  xi <- seq_len(nt)
  graphics::lines(xi, x, lwd=lwd, col=col, ...)
}

#' @rdname tapers-methods
#' @aliases points.tapers
#' @export
points.tapers <- function(x, pch="_", cex=1, ...){
  stopifnot(is.tapers(x))
  nt <- length(x)
  xi <- seq_len(nt)
  graphics::points(xi, x, pch=pch, cex=cex, ...)
}

#' @rdname tapers-methods
#' @aliases plot.tapers
#' @export
plot.tapers <- function(x, xi=NULL, color.pal=c("Blues","Spectral"), ylim=NULL, hv.lines=FALSE, log.y=FALSE, ...){
  stopifnot(is.tapers(x))
  nt <- length(x)
  if (is.null(xi)){
    xi <- seq_len(nt)
  }
  stopifnot(length(xi)==nt)
  mx <- max(x, na.rm=TRUE)
  pal <- match.arg(color.pal)
  # Color palette
  npal <- switch(pal, RdYlBu=11, Spectral=11, Blues=9)
  pal.col <- RColorBrewer::brewer.pal(npal, pal)
  # and turn into a function
  PALCOL <- grDevices::colorRampPalette(pal.col)
  cols <- PALCOL(mx)
  if (is.null(ylim)) ylim <- c(1, 1.1*mx)
  
  graphics::plot.default(xi, x,
                         ylab = "number of tapers",
                         xlab = "taper index",
                         ylim = ylim, 
                         yaxs = "i", xaxs = "i",
                         lwd = 1.8,
                         type = "h",
                         col = cols[x],
                         log = ifelse(log.y, "y", ""),
                         ...)
  
  graphics::lines(xi, x, col="black", lwd=0.7)
  
  if (hv.lines){
    # plot log2 multiples as horiz lines
    hl <- 2 ** (seq_len(ceiling(log2(mx))))
    graphics::abline(h=hl, lty=1, lwd=0.6, col="black")
    vl <- c(1, nt)
    graphics::abline(v=vl, lty=3, lwd=2, col="blue")
  }
  return(invisible(list(x=xi, k=x, line.colors=cols)))
}

###
###  Weighting methods
###

#' Calculate parabolic weighting factors.
#'
#' The resampled spectrum involves summing weighted tapers; this produces
#' the weighting factors. 
#' \code{\link{parabolic_weights_rcpp}} is the fastest implementation, used by
#' \code{\link{resample_fft_rcpp}}, but it takes only a single value.
#' \code{\link{parabolic_weights}} calls \code{\link{parabolic_weights_fast}} for vectors.
#'
#' If one has a \code{tapers} object, specify the \code{taper.index} to
#' produce a sequence of weights up to the value at that index; the user
#' is likely to never need to use this function though.
#'
#' Weighting factors, \eqn{W}{w}, are calculated as follows:
#' \deqn{
#'  W \equiv \frac{6 (n^2 - K^2)}{n (4 * n - 1) (n + 1)}
#' }{
#'  w = 6 (k^2 - (K-1)^2) / (k (4 * k - 1) (k + 1))
#' }
#' where \eqn{n}{k} is the total number of tapers, and 
#' \eqn{K}{K} is the integer sequence \eqn{[0,n-1]}{[0,K-1]}.
#' 
#' The sum of tapers should equal 1, within machine precision, when \eqn{n>0}{k>0}.
#'
#' @export
#' @author A.J. Barbour adapted the original algorithm (R.L. Parker), and authored the optimized versions.
#' @seealso \code{\link{resample_fft_rcpp}}, \code{\link{psdcore}}, \code{\link{riedsid}}
#'
#' @param tapvec \code{'tapers'} object; the number of tapers at each frequency
#' @param tap.index integer; the index of \code{tapvec} from which to produce a sequence of weights for
#' @param ntap integer; the number of tapers to provide weightings for.
#' 
#' @return A list with the number of tapers, indices of the taper sequence, and the weights \eqn{W_N}{w}.
#'
#' @example inst/Examples/rdex_parabolicweights.R
parabolic_weights <- function(tapvec, tap.index=1L) UseMethod("parabolic_weights")

#' @rdname parabolic_weights
#' @export
parabolic_weights.tapers <- function(tapvec, tap.index=1L){
  stopifnot(is.tapers(tapvec) | ((tap.index > 0L) & (tap.index <= length(tapvec))))
  kWeights <- parabolic_weights_fast(tapvec[as.integer(tap.index)])
  return(kWeights)
}

#' @rdname parabolic_weights
#' @export
parabolic_weights_fast <- function(ntap=1L) {
  # Must be long-integer, otherwise overflow for ntap > 1e3
  K <- as.double(ntap)
  kseq <- seq_len(ntap) - 1
  #
  ksq <- kseq * kseq # vector
  K2 <- K * K   	   # scalar
  K3 <- K2 * K 		   # scalar
  #
  # orig: w = (tapers^2 - (k-1).^2) * (1.5/(tapers*(tapers-0.25)*(tapers+1)));
  # or:       (tapers^2 - (k-1).^2) * 3/(2*K3 + K2*3/2 - K/2)
  # or:
  return(list(ntap=ntap, taper_seq = kseq + 1, taper_weights = (K2 - ksq) * 6 / ( 4 * K3  +  3 * K2  -  K )))
}

###
###  Constraint methods
###

#' @title Taper constraint methods
#'
#' @description
#' In the Riedel-Sidorenko recipe, the number of optimal tapers
#' at each frequency is strongly dependent on the first and
#' second derivatives of the spectrum. It is crucial to enforce
#' constraints on the number of actual tapers applied; this is
#' because the derivatives of "noisy" series can be bogus.
#' 
#' \code{\link{constrain_tapers}} refines the number of tapers at each frequency.
#' 
#' \code{\link{minspan}} sets bounds on the number of tapers at each frequency.
#' 
#' @details The method by which \code{\link{constrain_tapers}} refines tapers is 
#' set with the \code{constraint.method} argument:
#' 
#' \itemize{
#'   \item \code{'simple.slope.rcpp'} uses \code{\link{ctap_simple_rcpp}}
#'   \item \code{'simple.slope'} uses \code{\link{ctap_simple}}
#'   \item \code{'loess.smooth'} uses \code{\link{ctap_loess}}
#'   \item \code{'none'} returns unbounded tapers.
#' }
#' 
#' \code{\link{minspan}} bounds the number of tapers to within
#' the minimum of either the maximum number of tapers found in the object, 
#' or the half-length of the series, which is necessary because 
#' it would be nonsense to have more tapers than the length of the series. 
#' 
#' Details of the constraint methods:
#' 
#' \subsection{via first differencing (the default)}{
#' 
#' \code{\link{ctap_simple_rcpp}} is the preferred constraint method
#' (in previous versions \code{\link{ctap_simple}} was).
#' The algorithm uses first-differencing to modify the number
#' of tapers in the previous position.  Effectively, the constraint
#' is based on a causal, 1st-order Finite Impulse-response Filter (FIR) 
#' which makes the method sensitive to rapid changes in the number of tapers; 
#' naturally, smoother spectra tend to produce less fluctuation in taper numbers, 
#' which makes this well suited for adaptive processing. 
#'
#' This produces, generally, the most
#' stable results, meaning repeatedly running the constraint will not change values
#' other than on the first execution; the same cannot be said for the other
#' methods, which are also considerably more expensive to use.
#' 
#' }
#' 
#' \subsection{via LOESS smoothing}{
#' 
#' \code{\link{ctap_loess}} uses \code{\link{loess}} to smooth the taper vector; is
#' can be very slow thanks to quadratic scaling.
#' 
#' }
#'
#' @section Warning:
#'
#' \code{\link{ctap_loess}} results tend to be strongly dependent on
#' the tuning parameters given to \code{loess} (for obvious reasons); hence, 
#' some effort should be given to understand their effect, and/or re-tuning them if needed.
#'
#' @param tapvec \code{'tapers'} object; the number of tapers at each frequency
#' @param tapseq vector; positions or frequencies -- necessary for smoother methods
#' @param constraint.method  character; method to use for constraints on tapers numbers
#' @param verbose logical; should warnings and messages be given?
#' @param Kmin numeric; the minimum to set; default is 1
#' @param Kmax numeric; the maximum to set; default is the minimum of either (7/5 max value), or (1/2 series length)
#' @param ... optional arguments sent to the constrain function (e.g. \code{\link{ctap_simple}})
#' 
#' @return \code{\link{constrain_tapers}}: an object with class \code{'tapers'}; \code{\link{minspan}}: a vector
#' 
#' @author A.J. Barbour and R.L. Parker
#' 
#' @rdname tapers-constraints
#' @name tapers-constraints
#' 
#' @seealso \code{\link{riedsid}}, \code{\link{ctap_simple_rcpp}}, \code{\link{ctap_loess}}, \code{\link{tapers}}
#' @example inst/Examples/rdex_constraintapers.R
NULL

#' @rdname tapers-constraints
#' @export
constrain_tapers <- function(tapvec, ...) UseMethod("constrain_tapers")

#' @rdname tapers-constraints
#' @export
constrain_tapers.tapers <- function(tapvec, ...){
  as.tapers(constrain_tapers(as.vector(tapvec), ...), setspan=TRUE)
}

#' @rdname tapers-constraints
#' @export
constrain_tapers.default <- function(tapvec, tapseq=NULL,
                                     constraint.method=c("simple.slope.rcpp",
                                                         "simple.slope",
                                                         "loess.smooth",
                                                         "none"),
                                     verbose=TRUE, ...){
  # choose the appropriate method to apply taper constraints
  cmeth <- match.arg(constraint.method)
  tapvec.adj <- if (cmeth=="none"){
    if (verbose) warning("no taper optimization constraints applied")
    tapvec
  } else {
    if (verbose) message(sprintf("Constraining tapers with  ...  %s  ...  method", cmeth))
    if (cmeth == 'simple.slope.rcpp'){
      ctap_simple_rcpp(tapvec, ...)
    } else if (cmeth == 'simple.slope'){
      ctap_simple(tapvec, ...)
    } else if (cmeth == 'loess.smooth'){
      ctap_loess(tapvec, tapseq=tapseq, ...)
    } else {
      stop('no constraint function available')
    }
  }
  return(tapvec.adj)
}

#' @rdname tapers-constraints
#' @export
minspan <- function(tapvec, ...) UseMethod("minspan")

#' @rdname tapers-constraints
#' @export
minspan.tapers <- function(tapvec, ...){
  stopifnot(is.tapers(tapvec))
  tapvec.adj <- minspan(as.vector(tapvec))
  return(as.tapers(tapvec.adj))
}

#' @rdname tapers-constraints
#' @export
minspan.default <- function(tapvec, Kmin=NULL, Kmax=NULL, ...){
  tapvec <- as.integer( tapvec )
  if (is.null(Kmin)) Kmin <- 1
  if (is.null(Kmax)){
    Kmax.upper.a <- floor(7 * max(tapvec, na.rm=TRUE) / 5)
    Kmax.upper.b <- floor(length(tapvec) / 2)
    Kmax <- min(c(Kmax.upper.a, Kmax.upper.b), na.rm=TRUE)
  }
  stopifnot(Kmin <= Kmax)
  tapvec <- as.integer( pmax(pmin(tapvec, Kmax, na.rm=TRUE), Kmin, na.rm=TRUE) )
  stopifnot(all(tapvec >= 0))
  return(tapvec)
}

##
## Individual methods:
##

#' @title Taper constraints using simple derivatives
#' @details
#' \code{\link{ctap_simple}} is the original version ported to c, and
#' \code{\link{ctap_simple_rcpp}} is the recommended version to use.
#' @export
#' @param tapvec integer; the number of tapers at each frequency (can be a vector)
#' @param maxslope integer; constrain based on this maximum first difference
#' @param ... additional arguments
#' @seealso \code{\link{constrain_tapers}}, \code{\link{ctap_loess}}
#' @examples
#' 
#' # generate some random taper series and constrain them based on slopes
#' set.seed(1237)
#' n <- 11
#' x <- seq_len(n)
#' xn <- round(runif(n,1,n))
#' 
#' xnf <- ctap_simple_rcpp(xn, 0) # flattens out
#' xnc <- ctap_simple_rcpp(xn, 1) # no change, already only slopes = 1
#' try(all.equal(xnc, xn))
#' xnc2 <- ctap_simple_rcpp(xn, 2) # slopes = 2 only
#'
#' plot(xn, type='b', pch=16, ylim=c(0,12))
#' grid()
#' abline(a=0,b=1, col='red', lty=3); abline(a=0,b=2, col='blue', lty=3)
#' lines(xnf, type='b', col='green')
#' lines(xnc, type='b', col='red')
#' lines(xnc2, type='b', col='blue')
#' lines(0.2+as.vector(psd::ctap_simple(psd::as.tapers(xn))), type='b', pch=".", col='salmon')
#'
#' # compare simple and rcpp implementations
#' kcr <- ctap_simple_rcpp(xn, 2)
#' kcs <- ctap_simple(xn, 2)
#' rbind(kcs, kcr)
#' try(all.equal(kcr, kcs))
#'
#' # more examples:
ctap_simple_rcpp <- function(tapvec, ...) UseMethod("ctap_simple_rcpp")

#' @rdname ctap_simple_rcpp
#' @export
ctap_simple_rcpp.tapers <- function(tapvec, ...){
  # c++ code used for speed up of forward+backward operations
  tapvec <- as.integer(tapvec)
  tapvec.adj <- ctap_simple_rcpp.default(tapvec, ...)
  return(as.tapers(tapvec.adj))
}

#' @rdname ctap_simple_rcpp
#' @export
ctap_simple_rcpp.default <- function(tapvec, maxslope = 1L, ...) {
  tapvec <- as.integer(tapvec)
  maxslope <- as.integer(maxslope)
  tapvec.adj <- rcpp_ctap_simple(tapvec, maxslope)
  return(tapvec.adj)
}

#' @rdname ctap_simple_rcpp
#' @export
ctap_simple <- function(tapvec, ...) UseMethod("ctap_simple")

#' @rdname ctap_simple_rcpp
#' @export
ctap_simple.tapers <- function(tapvec, ...){
  stopifnot(is.tapers(tapvec))
  tapvec.adj <- ctap_simple.default(as.vector(tapvec), ...)
  return(as.tapers(tapvec.adj))
}

#' @rdname ctap_simple_rcpp
#' @export
ctap_simple.default <- function(tapvec, maxslope=1L, ...){
  # current code requires real
  tapvec <- as.numeric(tapvec)
  maxslope <- as.numeric(maxslope)
  # c code used for speed up of forward+backward operations
  tapvec.adj <- as.integer(.Call("rlp_constrain_tapers", tapvec, maxslope, PACKAGE="psd"))
  return(tapvec.adj)
}

#' @title Taper constraints using loess smoothing
#' @rdname ctap_loess
#' @export
#' @param tapvec integer; the number of tapers at each frequency (can be a vector)
#' @param tapseq vector; positions to evaluate derivatives (unused here, but necessary for smoother methods)
#' @param loess.span  scalar; the span used in \code{loess}
#' @param loess.degree  scalar; the polynomial degree
#' @param verbose logical; should warnings and messages be given?
#' @param ... additional arguments
#' @seealso \code{\link{constrain_tapers}}, \code{\link{ctap_simple_rcpp}}
ctap_loess <- function(tapvec, ...) UseMethod("ctap_loess")

#' @rdname ctap_loess
#' @export
ctap_loess.tapers <- function(tapvec, ...){
  stopifnot(is.tapers(tapvec))
  tapvec.adj <- ctap_loess(as.vector(tapvec), ...)
  return(as.tapers(tapvec.adj))
}


#' @rdname ctap_loess
#' @export
ctap_loess.default <- function(tapvec, tapseq=NULL, loess.span=.3, loess.degree=1, verbose=TRUE, ...){
  # having an appropriate x-sequence is absolutely critical to obtaining useful results
  if (is.null(tapseq)){
    tapseq <- seq_along(tapvec)
    if (verbose) warning("Generated a position sequence; results may be bogus.")
  }
  lt <- length(tapvec)
  if (verbose & lt > 1e4) warning("Loess-method has quadratic memory scaling (1e3 pt -> 10 Mb)...")
  trc <- ifelse(lt >= 1e3, "approximate", "exact")
  loe <- stats::loess(y ~ x, 
                      data.frame(x=tapseq, y=as.numeric(tapvec)), 
                      span=loess.span, degree=loess.degree,
                      control = loess.control(trace.hat = trc))
  tapvec.adj <- as.integer(stats::predict(loe))
  return(tapvec.adj)
}
