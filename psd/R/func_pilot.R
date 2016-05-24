#' Calculate inital power spectral density estimates
#'
#' @description
#' This PSD is used as the starting point -- the pilot spectrum -- for
#' the adaptive estimation routine.
#'
#' @details
#' A fixed number
#' of tapers is applied across all frequencies using \code{\link{psdcore}}, and
#' subsequent taper-refinements are based on the spectral derivatives
#' of this spectrum; hence, changes in the number of tapers can affect
#' how many adaptive stages may be needed (though there are no formal convergence
#' criteria to speak of).
#'
#' The taper series of the returned spectrum is constrained using
#' \code{as.tapers(..., minspan=TRUE)}.
#'
#' The default behaviour (\code{remove.AR <= 0}) is to remove the standard linear 
#' model \eqn{[f(x) = \alpha x + \beta]} from the data; however,
#' the user can model the effect of an autoregressive process by specifiying
#' \code{remove.AR}.
#'
#' @section Removing an AR effect from the spectrum:
#' If \code{remove.AR > 0} the argument is used as \code{AR.max} in 
#' \code{\link{prewhiten}}, from which an AR-response spectrum is calculated using
#' the best fitting model.
#'
#' If the value of \code{remove.AR} is too low the spectrum 
#' could become distorted,
#' so use with care.
#' \emph{Note, however, that the 
#' value of \code{remove.AR} will be restricted to within the 
#' range \eqn{[1,100]}.}
#' If the AR order is much larger than this, it's unclear how \code{\link{prewhiten}}
#' will perform and whether the AR model is appropriate.
#'
#' \emph{Note that this function does not produce a parametric spectrum estimation; rather,
#' it will return the amplitude response of the best-fitting AR model as \code{\link[stats]{spec.ar}}
#' would. \strong{Interpret these results with caution, as an AR response spectrum
#' can be misleading.}}
#'
#' @name pilot_spec
#' @aliases pilot_spectrum spec.pilot
#' @export
#' @author A.J. Barbour
#' @seealso \code{\link{psdcore}}, \code{\link{prewhiten}}, \code{\link[stats]{spec.ar}}
#'
#' @param x  vector; the data series to find a pilot spectrum for
#' @param x.frequency  scalar; the sampling frequency (e.g. Hz) of the series
#' @param ntap  scalar; the number of tapers to apply during spectrum estimation
#' @param remove.AR  scalar; the max AR model to be removed from the data.
#' @param plot  logical; should a plot be created?
#' @param verbose  logical; should messages be given?
#' @param ...  additional parameters passed to \code{\link{psdcore}}
#' @return An object with class 'spec', invisibly, and \code{"pilot_psd"} in the working environment.
#'
#' @example inst/Examples/rdex_pilotspec.R
pilot_spec <- function(x, ...) UseMethod("pilot_spec")

#' @rdname pilot_spec
#' @aliases pilot_spec.ts
#' @export
pilot_spec.ts <- function(x, ...){
  stopifnot(is.ts(x))
  frq <- frequency(x)
  pilot_spec.default(as.vector(x), x.frequency=frq, ...)  
}

#' @rdname pilot_spec
#' @aliases pilot_spec.default
#' @export
pilot_spec.default <- function(x, x.frequency=NULL, ntap=NULL, remove.AR=NULL, plot=FALSE, verbose=FALSE, ...){
  
  if (is.null(ntap)) ntap <- 7
  if (is.null(remove.AR)) remove.AR <- 0
  if (is.null(x.frequency)) x.frequency <- 1
  stopifnot(length(ntap)==1)
  stopifnot(length(remove.AR)==1)
  stopifnot(length(x.frequency)==1)

  # AR spectrum or no?
  REMAR <- ifelse(remove.AR > 0, TRUE, FALSE)
  
  #restrict maximum ar orders to within [1,100]
  if (REMAR) remove.AR <- max(1, min(100, abs(remove.AR)))
  
  xprew <- prewhiten(x, x.fsamp=x.frequency, AR.max=remove.AR, detrend=TRUE, impute=TRUE, plot=FALSE, verbose=verbose)
  
  ## Remove and AR model
  if (REMAR){
    # AR fit
    ordAR <- xprew[['ardfit']][['order']]
    if (ordAR==0){
      warning("AR(0) was the highest model found!\n\t\tConsider fitting a linear model instead ( remove.AR = 0 ).")
    } else {
      if (verbose) message(sprintf("removed AR(%s) effects from the spectrum", ordAR))
    }
    # calculate PSD of the AR fit
    xar <- xprew[['prew_ar']] # ts object
    Pspec_ar <- psdcore(xar, ntaper=ntap, AR=TRUE, preproc=FALSE, refresh=TRUE, verbose=FALSE)
    arvar <- var(Pspec_ar[['spec']])
    mARs <- mean(Pspec_ar[['spec']])
    Pspec_ar[['spec']] <- Pspec_ar[['spec']] / mARs
  }
  
  ## Calculate spectrum of non-AR (lm) model:
  xlm <- xprew[['prew_lm']] # ts object
  Pspec <- psdcore(xlm, ntaper=ntap, AR=FALSE, preproc=FALSE, refresh=TRUE, verbose=FALSE)
  Ptap <- Pspec[['taper']]
  num_tap <- length(Ptap)
  num_frq <- Pspec[['numfreq']]
  
  ## generate a series, if necessary
  if (num_tap < num_frq) Ptap <- rep.int(Ptap[1], num_frq)
  
  ## return tapers object
  Pspec[['taper']] <- as.tapers(Ptap, setspan=TRUE)
  
  ## remove the spectrum of the AR process
  if (REMAR){
    stopifnot(exists("Pspec") & exists("Pspec_ar"))
    if (verbose) message(sprintf("Removing AR(%s) effects from spectrum", ordAR))
    # reup the spectrum
    Ospec <- psd_envAssignGet("pre_AR_psd", Pspec)
    psd_envAssign("AR_psd", Pspec_ar)
    p.lin <- Pspec[['spec']]
    p.ar <- Pspec_ar[['spec']]
    Pspec[['spec']] <-  p.lin / p.ar
  }
  
  if (plot){
    ttl <- "Pilot spectrum estimation"
    llog <- 'dB'
    try({
      if (REMAR){
        if (verbose) message('Plotting,', tolower(ttl))
        par(las=1)
        plot(Ospec, log=llog, col="red", main=ttl)
        mtext(sprintf("(with AR(%s) response)", ordAR), line=0.4)
        # rescale
        Pspec_ar[['spec']] <- Pspec_ar[['spec']] * mARs
        plot(Pspec_ar, log=llog, col="blue", add=TRUE)
        plot(Pspec, log=llog, add=TRUE, lwd=2)
        legend("bottomleft", 
               c("original PSD",
                 sprintf("AR-innovations PSD\n(mean %.01f +- %.01f dB)", dB(mARs), dB(sqrt(arvar))/4),
                 "AR-filter response"), 
               lwd=2, col=c("red","blue","black"))
      } else {
        plot(Pspec, log=llog, main=ttl)
      }
    })
  }
  return(invisible(psd_envAssignGet("pilot_psd", Pspec)))
}

