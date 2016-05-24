#' Adaptive sine multitaper power spectral density estimation
#' 
#' @description
#' This is the primary function to be used in this package: it returns
#' power spectral density estimates of a univariate timeseries, with
#' an optimal number of tapers at each frequency based on iterative
#' reweighted spectral derivatives.
#'
#' @details
#' See the \strong{Adaptive estimation} section in the description of
#' the \code{\link{psd-package}} for details regarding adaptive estimation.
#' 
#' \code{\link{pspectrum_basic}} is a simplified implementation used mainly for
#' testing.
#'
#' @name pspectrum
#' @export
#' @author A.J. Barbour adapted original by R.L. Parker
#' @seealso \code{\link{psdcore}}, \code{\link{pilot_spec}}, \code{\link{riedsid}}, \code{\link{prewhiten}}
#' 
#' @param x vector; series to find PSD estimates for
#' @param x.frqsamp scalar; the sampling rate (e.g. Hz) of the series \code{x}; equivalent to \code{\link{frequency}}.
#' @param ntap.init scalar; the number of sine tapers to use in the pilot spectrum estimation; if \code{NULL} then the
#' default in \code{\link{pilot_spec}} is used.
#' @param niter scalar; the number of adaptive iterations to execute after the pilot spectrum is estimated.
#' @param AR logical; should the effects of an AR model be removed from the pilot spectrum?
#' @param Nyquist.normalize  logical; should the units be returned on Hz, rather than Nyquist?
#' @param verbose logical; Should messages be given?
#' @param no.history logical; Should the adaptive history \emph{not} be saved?
#' @param plot logical; Should the results be plotted?
#' @param ... Optional parameters passed to \code{\link{riedsid}}
#' @param stage integer; the current adaptive stage (0 is pilot)
#' @param dvar numeric; the spectral variance; see also \code{\link{vardiff}} etc
#' @return Object with class 'spec', invisibly. It also assigns the object to
#' \code{"final_psd"} in the working environment.
#'
#' @example inst/Examples/rdex_pspectrum.R
#' 
pspectrum <- function(x, ...) UseMethod("pspectrum")

#' @rdname pspectrum
#' @aliases pspectrum.ts
#' @export
pspectrum.ts <- function(x, ...){
  frq <- frequency(x)
  pspectrum(as.vector(x), x.frqsamp=frq, ...)  
}

#' @rdname pspectrum
#' @aliases pspectrum.spec
#' @export
pspectrum.spec <- function(x, ...){
  cant <- "cannot adapt  pspectrum  results without an fft in the psd environment. see ?pspectrum"
  if (inherits(x, "amt")){
    name <- getOption("psd.ops")[['names']]
    fft <- psd_envGet(name[['fft']])
    if (is.null(fft)){
      stop(cant)
    } else {
      warning('updating  pspectrum  results is not (yet) supported') 
      .NotYetImplemented()
    }
  } else {
    .NotYetImplemented()
  }
}

#' @rdname pspectrum
#' @aliases pspectrum.default
#' @export
pspectrum.default <- function(x, x.frqsamp=1, ntap.init=NULL, niter=5, AR=FALSE, Nyquist.normalize=TRUE, verbose=TRUE, no.history=FALSE, plot=FALSE, ...){
  
  stopifnot(length(x)>1)
  
  # plotting and iterations
  if (is.null(niter)) stopifnot(niter>=0)
  plotpsd_ <- FALSE
  # iteration stages (0 is pilot)
  iter_stages <- 0:niter
  
  # retain history
  save_hist <- ifelse((niter < 10) & !no.history, TRUE, FALSE)
  
  # AR switch
  ordAR <- ifelse(AR, 100, 0)

  for (stage in iter_stages){
    
    if (stage==0){
      if (verbose) adapt_message(stage)
      if (niter==0 & plot) plotpsd_ <- TRUE
      # --- setup the environment ---
      psd_envRefresh(verbose=verbose)

      # --- pilot spec ---
      # ** normalization is here:
      Pspec <- pilot_spec(x, x.frequency=x.frqsamp, ntap=ntap.init, 
                          remove.AR=ordAR, verbose=verbose, plot=plotpsd_)
      kopt <- Pspec[['taper']]
      
      # ensure series is in the environment
      psd_envAssign("original_pspectrum_series", x)
      
      # starting spec variance
      dvar.o <- varddiff(Pspec)
      
      # --- history ---
      if (save_hist){
        new_adapt_history(niter)
        update_adapt_history(Pspec, stage)
      }
      
      x <- 0 # to prevent passing orig data back/forth
      
    } else {
      
      # enforce silence in the subfunctions once the adapting gets going
      rverb <- ifelse(stage > 0, FALSE, TRUE)
      
      ## calculate optimal tapers
      kopt <- riedsid2(Pspec, verbose=rverb, ...)
      
      # get data back for plotting, etc.
      if (stage==niter){
        x <- psd_envGet("original_pspectrum_series")
        if (plot){
          plotpsd_ <- TRUE
        }
      }
  
      # update spectrum with new tapers
      # TODO: here's why preproc flags are wrong...
      Pspec <- psdcore(X.d=x, X.frq=x.frqsamp, ntaper=kopt, 
                       preproc=FALSE, plot=plotpsd_, verbose=rverb) 
      
      # show spectral variance reduction
      if (verbose) adapt_message(stage, varddiff(Pspec)/dvar.o)
      
      ## update history
      if (save_hist) update_adapt_history(Pspec, stage)
      
    }
  }
  if (Nyquist.normalize) Pspec <- normalize(Pspec, x.frqsamp, verbose=verbose)
  return(invisible(psd_envAssignGet("final_psd", Pspec)))
}

#' @rdname pspectrum
#' @export
pspectrum_basic <- function(x, ntap.init=7, niter=5, verbose=TRUE, ...){
  
  if (verbose) adapt_message(0)
  P <- psdcore(x, ntaper=ntap.init, preproc = FALSE, refresh=TRUE)
  ko <- P[['taper']]
  nf <- length(ko)
  
  # Iterate on optimal tapers, and resample spectrum
  if (verbose & niter > 0) message("Iterative refinement of spectrum (", niter, " iterations)")
  for (iter in seq_len(niter)){
    if (verbose) adapt_message(iter)
    # find optimal tapers
    ko <- riedsid2(P, ko, verbose=FALSE)
    # update spectrum
    P  <- psdcore(x, ntaper=ko, preproc = FALSE)
  }
  return(P)
}

#' @rdname pspectrum
#' @export
adapt_message <- function(stage, dvar=NULL){
  stopifnot(stage >= 0)
  stage <- if (stage == 0){
    paste(stage, "est. (pilot)")
  } else {
    if (!is.null(dvar)){
      paste(stage, sprintf("est. (Ave. S.V.R. %.01f dB)", dB(dvar)))
    } else {
      stage
    }
  }
  message(sprintf("Stage  %s ", stage))
}

