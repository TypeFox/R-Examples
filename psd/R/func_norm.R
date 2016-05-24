#' Normalization of power spectral density estimates.
#'
#' Normalize power spectral densities from
#' various estimators into single-sided spectra.
#'
#' Normalizations commonly encountered for power spectra 
#' depend on it's assumed sidedness: whether the
#' spectrum is either single- or double-sided.
#' The normalizations performed here enforce single-sidedness, and correct
#' as necessary.
#'
#' Frequencies are assumed to be based on the Nyquist frequency (half the 
#' sampling rate).  For example: If a series \eqn{X} has sampling frequency \eqn{F_S},
#' then the PSD frequencies will span \eqn{[0,F_S/2]}.
#'
#' For amplitudes, improper normalization can can introduce errant factors
#' of either 1/2 or \eqn{F_S} into the estimates, depending on the assumed sidedness.  
#' These factors can be accounted for with the \code{src}
#' argument, which defaults to normalizing a double-sided spectrum.
#' 
#' @section Spectrum sidedness and the \code{src} argument:
#' \subsection{\code{"double.sided"} or \code{"spectrum"}}{
#'
#' These spectra assume frequency range of \eqn{[-F_S/2,F_S/2]}, and so are normalized
#' by scaling by a factor of two upwards.
#' Some estimators producing double-sided spectra: 
#' \itemize{
#' \item{\code{stats::spectrum}}{}
#' \item{\code{RSEIS::mtapspec}}{}
#' }
#' }
#'
#' \subsection{\code{"single.sided"} or \code{"psd"}}{
#' As mentioned before, 
#' these spectra assume frequency range of \eqn{[0,F_S/2]} and
#' are scaled only by the inverse of the sampling rate.
#' Some estimators producing single-sided spectra: 
#' \itemize{
#' \item{\code{\link{psdcore}}}{}
#' }
#' }
#'
#' @name psd-normalization
#' @rdname psd-normalization
#' @aliases normalization
#' @author A.J. Barbour
#'
#' @param Spec spectrum to normalize
#' @param Fsamp sampling frequency
#' @param src character string; the source of the spectrum estimator
#' @param verbose logical; should messages be given?
#' @param ... (unused) additional parameters
#' @return An object with its spectral values normalized accordingly.
#'
#' @seealso \code{\link{psdcore}}, \code{\link{spectral_properties}}
#' @example inst/Examples/rdex_normalization.R
NULL
 
#' @rdname psd-normalization
#' @aliases normalize
#' @export
normalize <- function(Spec, ...) UseMethod("normalize")

#' @rdname psd-normalization
#' @export
normalize.list <- function(Spec, ...){
  stopifnot(exists("freq", where=Spec) & exists("spec", where=Spec))
  class(Spec) <- 'spec'
  Spec <- normalize(Spec, ...)
  class(Spec) <- 'list'
  return(Spec)
}

#' @rdname psd-normalization
#' @aliases normalize.spec
#' @export
normalize.spec <- function(Spec, Fsamp=1, src=c("spectrum","double.sided","psd","single.sided"), verbose=TRUE, ...){
  
  stopifnot(is.spec(Spec))
  
  Fsamp <- if (Fsamp > 0){
    # value represents sampling frequency
    Fsamp
  } else if (Fsamp < 0){
    # value is sampling interval
    abs(1/Fsamp)
  } else {
    stop("bad sampling information")
  }
  
  # Find out which type to apply
  src <- match.arg(src)
  is.single.sided <- switch(toupper(src), SPECTRUM=FALSE, DOUBLE.SIDED=FALSE, PSD=TRUE, SINGLE.SIDED=TRUE)
  ptyp <- if (is.single.sided){
    # spectrum is from psd, and is single-sided
    Spec[['spec']] <- Spec[['spec']] / Fsamp
    "single"
  } else {
    # spectrum is from spectrum or others, double sided
    Spec[['spec']] <- Spec[['spec']] * 2
    "double"
  }
  
  if (verbose) message(sprintf("Normalized  %s-sided  psd estimates ( %s ) for sampling-freq.  %s", ptyp, src, Fsamp))
  return(invisible(Spec))
}

#' @rdname psd-normalization
#' @aliases normalize.amt
#' @export
normalize.amt <- function(Spec, ...){
  Spec <- normalize.spec(Spec, src='psd', ...)
  return(Spec)
}

