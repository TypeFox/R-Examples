#' Constrained, optimal tapers using the Riedel & Sidorenko--Parker method
#' 
#' Estimates the optimal number of tapers at each frequency of
#' given PSD, using a modified Riedel-Sidorenko MSE recipe (RS-RLP).
#' 
#' @details
#' The optimization is as follows. First, weighted derivatives of the 
#' input PSD are computed.
#' Using those derivates the optimal number of tapers is found through the 
#' RS-RLP formulation.
#' Constraints are then placed on the practicable number of tapers.
#' 
#' \code{\link{riedsid2}} is a new implementation which does not allow 
#' for multiple constraint methods; this is the preferred function to use.
#'
#' \subsection{Taper constraints}{
#' The parameter \code{c.method} provides an option to change the method
#' of taper constraints.  A description of each may be found in 
#' the documentation for \code{\link{constrain_tapers}}.
#'
#' Once can use \code{constrained=FALSE} to turn off all taper constraints; this
#' could lead to strange behavior though.
#' }
#'
#' \subsection{Spectral derivatives}{
#' The parameter \code{Deriv.method} determines which method is used
#' to estimate derivatives.
#' \itemize{
#' \item{\code{"local_qls"}}{ (\strong{default}) uses quadratic weighting and
#' local least-squares estimation; this can be slower than \code{"spg"}.}
#' \item{\code{"spg"}}{ uses \code{\link{splineGrad}}; then, additional arguments
#' may be passed to control the smoothness of the derivatives
#' (e.g \code{spar} in \code{smooth.spline}).}
#' }
#' }
#'
#' @section Warning:
#' The \code{"spg"} can become numerically unstable, and it's not clear when it will
#' be the preferred over the \code{"local_qls"} method, other than for efficiency's sake.
#'
#' @export
#' @author A.J. Barbour adapted original by R.L. Parker
#' 
#' @param PSD vector or class \code{'amt'} or \code{'spec'}; the spectral values used to optimize taper numbers
#' @param ntaper scalar or vector; number of tapers to apply optimization
#' @param tapseq vector; representing positions or frequencies (same length as \code{PSD})
#' @param Deriv.method character; choice of gradient estimation method 
#' @param constrained logical; apply constraints with \code{\link{constrain_tapers}}; \code{FALSE} turns off constraints
#' @param c.method string; constraint method to use with \code{\link{constrain_tapers}}, only if \code{constrained=TRUE}
#' @param verbose logical; should messages be printed?
#' @param ... optional argments passed to \code{\link{constrain_tapers}}
#' @return Object with class \code{'tapers'}
#' 
#' @seealso \code{\link{constrain_tapers}}, \code{\link{resample_fft_rcpp}}, \code{\link{psdcore}}, \code{\link{pspectrum}}
#' @example inst/Examples/rdex_riedsid.R
riedsid <- function(PSD, ...) UseMethod("riedsid")

#' @rdname riedsid
#' @export
riedsid.spec <- function(PSD, ...){
  stopifnot(is.spec(PSD))
  Pspec <- PSD[['spec']]
  Tapseq <- PSD[['freq']]
  ntaper <- if (is.amt(PSD)){
      PSD[['taper']]
  } else {
      rep.int(1L, length(Pspec))
  }
  riedsid.default(PSD=Pspec, ntaper=ntaper, tapseq=Tapseq, ...)
}


#' @rdname riedsid
#' @export
riedsid.default <- function(PSD, ntaper = 1L, 
                            tapseq=NULL, 
                            Deriv.method=c("local_qls","spg"),
                            constrained=TRUE, c.method=NULL,
                            verbose=TRUE, ...) {
  ## spectral values
  PSD <- as.vector(PSD)
  # num freqs
  nf <- psd_envAssignGet("num_freqs", length(PSD))
  # prelims
  #  A small number to protect against zeros
  eps <- 1e-78 
  # vectorize initial estimate
  Zeros <- zeros(nf)
  nt <- length(ntaper)
  ntap <- if (nt==1){
    ntaper + Zeros
  } else {
    ntaper
  }
  if (!is.tapers(ntap)) ntap <- as.tapers(ntap)
  
  # Set the number of tapers to within the range: 1/2 nf, 7/5 ntap
  # rowMins produces a rowvec of rowwise minimums; convert to colvec
  nspan <- minspan(ntap, nf)
  
  # The spectral gradients should be in log-space, so
  # create a log spec, and pad to handle begnning and end values
  nadd <- 1 + max(nspan)
  Y <- c(PSD[nadd:2], PSD, PSD[(nf-1):(nf-nadd)])
  Y[Y <= 0] <- eps
  lY <- log(Y)
  dY <- d2Y <- Zeros
  #
  kseq <- if (is.null(tapseq) | (length(tapseq) != nf)){
    seq.int(from=0, to=1/2, length.out=length(PSD))
  } else {
    tapseq
  }
  #
  # Smooth spectral derivatives
  #
  lsmeth <- switch(match.arg(Deriv.method), local_qls=TRUE, spg=FALSE)
  stopifnot(exists("lsmeth"))
  rss <- if (lsmeth){
    # spectral derivatives the preferred way
    DFUN <- function(j, 
                     j1=j-nspan[j]+nadd-1, 
                     j2=j+nspan[j]+nadd-1, 
                     jr=j1:j2, 
                     logY=lY[jr], 
                     dEps=eps){
      u <- jr - (j1 + j2)/2 # rowvec 
      u2 <- u*u             # rowvec
      L <- j2-j1+1          # constant
      L2 <- L*L             # constant
      LL2 <- L*L2           # constant
      LL2L <- LL2 - L       # constant
      #
      CC <- 12
      uzero <- (L2 - 1)/CC  # constant
      #
      # first deriv
      dY <- u %*% logY * CC / LL2L
      # second deriv
      d2Y <- (u2 - uzero) %*% logY * 360 / LL2L / (L2 - 4)
      return(c(fdY2=dY*dY, fd2Y=d2Y, fdEps=dEps))
    }
    DX <- seq_len(nf)
    ##
    ## Bottleneck:
    RSS <- vapply(X=DX, FUN=DFUN, FUN.VALUE=c(1,1,1))
    ##
    attr(RSS, which="lsderiv") <- lsmeth
    RSS <- psd_envAssignGet("spectral_derivatives.ls", RSS)
    msg <- "local quadratic regression"
    # sums:
    #[ ,1] fdY2
    #[ ,2] fd2Y
    #[ ,3] fdEps
    abs(colSums(RSS))
  } else {
    RSS <- splineGrad(dseq=log(0.5+kseq), #seq.int(0,.5,length.out=length(PSD))), 
                      dsig=log(PSD),
                      plot.derivs=FALSE, ...) #, spar=1)
    attr(RSS, which="lsderiv") <- lsmeth
    RSS <- psd_envAssignGet("spectral_derivatives", RSS) 
    #returns log
    RSS[,2:4] <- exp(RSS[,2:4])
    msg <- "weighted cubic spline"
    abs(eps + RSS[,4] + RSS[,3]**2)
  }
  if (verbose) message(sprintf("Using spectral derivatives from  %s", msg))
  #
  #(480)^0.2*abs(PSD/d2psd)^0.4
  # Original form:  kopt = 3.428*abs(PSD ./ d2psd).^0.4;
  # kopt = round( 3.428 ./ abs(eps + d2Y + dY.^2).^0.4 );
  #
  sc <- ifelse(TRUE, 473.3736, 480)
  kopt <- (sc ** 0.2)/(rss ** 0.4)
  #
  # Constrain tapers
  kopt <- if (constrained){
    constrain_tapers(tapvec = kopt, tapseq=kseq, constraint.method=c.method, verbose=verbose, ...)
  } else {
    as.tapers(kopt)
  }
  #
  return(kopt)
} 

#' @rdname riedsid
#' @export
riedsid2 <- function(PSD, ...) UseMethod("riedsid2")

#' @rdname riedsid
#' @export
riedsid2.spec <- function(PSD, ...){
  stopifnot(is.spec(PSD))
  pspec <- PSD[['spec']]
  freqs <- PSD[['freq']]
  ntaper <- if (is.amt(PSD)){
    PSD[['taper']]
  } else {
    rep.int(1L, length(pspec))
  }
  riedsid2.default(pspec, ntaper, ...)
}

#' @rdname riedsid
#' @export
riedsid2.default <- function(PSD, ntaper=1L, constrained=TRUE, verbose=TRUE, ...){
  
  PSD <- as.vector(PSD)
  ntaper <- as.vector(ntaper)
  
  #  A small number to protect against zeros
  eps <- 1e-78
  nf <- length(PSD)
  nt <- length(ntaper)
  if (nt == 1) ntaper <- rep(ntaper, nf)
  # some constraints
  nspan <- ceiling( pmin( nf/2, 7*ntaper/5 ) )
  nadd <- 1 + max(nspan)
  # Create log psd, and pad to handle beginning and end values
  ist <- nadd:2
  iend <- (nf - 1):(nf - nadd)
  S <- as.numeric(c(PSD[ist], PSD, PSD[iend])) + eps
  Y <- log(S)
  DFUN <- function(j){
    j1 <- j - nspan[j] + nadd - 1
    j2 <- j + nspan[j] + nadd - 1
    jseq <- j1:j2
    u <- jseq - (j1 + j2)/2
    L <- j2 - j1 + 1
    CC <- 12
    #
    uzero <- (L^2 - 1)/CC
    #
    # first deriv
    dY <- u  %*%  Y[jseq] * CC / (L*(L*L - 1))
    # second deriv
    d2Y <- (u*u - uzero)  %*%  Y[jseq] * 360 / (L*(L^2 - 1)*(L^2 - 4))
    #
    return(c(eps=eps, d2Y=d2Y, dYsq=dY*dY))
  }
  # Calculate derivatives
  yders <- vapply(X=seq_len(nf), FUN=DFUN, FUN.VALUE=c(1,1,1))
  # and optimal tapers
  sc <- ifelse(TRUE, 473.3736, 480)
  kopt <- round( sc**0.2 / abs(colSums(yders))**0.4 )
  
  kopt <- if (constrained){
    constrain_tapers(tapvec = kopt, verbose = verbose, ...)
  } else {
    as.tapers(kopt)
  }
  return(kopt)
}