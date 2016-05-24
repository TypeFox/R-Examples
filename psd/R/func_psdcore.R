#' Multitaper power spectral density estimates of a series
#'
#' Compute power spectral density (PSD) estimates
#' for the input series using sine multitapers.
#' This is used by \code{\link{pspectrum}} for the adaptive
#' estimation procedure.
#'  
#' @details
#' \subsection{Tapering}{
#' The parameter \code{ntaper} specifies the number of sine tapers to be used 
#' at each frequency: equal tapers at each frequency for a scalar; 
#' otherwise, use \code{ntaper[j]} sine tapers at \code{frequency[j]}.
#' }
#'
#' \subsection{Truncation}{
#' The series, with length \eqn{N}, is necessarily truncated so that \eqn{1+N/2} evenly 
#' spaced frequencies are returned. This truncation makes the series length ``highly composite",
#' which the discrete Fourier transform (DFT) is most efficient.
#' The "fftw" vignette (accessed with \code{vignette("fftw",package="psd")}) shows
#' how the performance of a DFT can be affected by series length.
#' }
#'
#' \subsection{Decimation}{
#' 	No longer supported. Setting \code{ndecimate} will not affect the results
#' }
#'
#' \subsection{Sampling}{
#'  If \code{X.frq} is NULL, the value is assumed to be 1, unless \code{X.d} is a  \code{'ts'} object.
#'  If \code{X.frq > 0} it's assumed the value represents \emph{frequency} (e.g. Hz).
#'  If \code{X.frq < 0} it's assumed the value represents \emph{interval} (e.g. seconds).
#' }
#'
#' @param X.d  the series to estimate a spectrum for 
#' @param X.frq  scalar; the sampling information (see section Sampling)
#' @param ntaper  scalar, vector, or \code{\link{tapers}}; the number of sine tapers to apply at each frequency
#' @param preproc  logical; should \code{X.d} have a linear trend removed?
#' @param na.action  function to deal with \code{NA} values
#' @param plot  logical; should the estimates be shown compared to the \code{\link[stats]{spectrum}}-based estimates?
#' Note that this will add some computation time, since the cosine-tapered periodogram is calculated inside
#' \code{\link{pgram_compare}}.
#' @param refresh  logical; ensure a free environment prior to execution
#' @param verbose logical; should warnings and messages be given?
#' @param ndecimate  now ignored
#' @param ... additional parameters
#' 
#' @return An on object of class \code{'amt','spec'}, which has a structure similar to a regular \code{'spec'} object, 
#' but with a few additional fields, invisibly.
#'
#' @name psdcore
#' @export
#' @author A.J. Barbour; original algorithm by R.L. Parker.
#' @seealso \code{\link{pspectrum}}, \code{\link{riedsid}}, \code{\link{parabolic_weights}}, \code{\link{pgram_compare}}
#'
#' @example inst/Examples/rdex_psdcore.R
psdcore <- function(X.d, ...) UseMethod("psdcore")

#' @rdname psdcore
#' @aliases psdcore.ts
#' @export
psdcore.ts <- function(X.d, ...){
  frq <- stats::frequency(X.d)
  psdcore(as.vector(X.d), X.frq = frq, ...)
}

#' @rdname psdcore
#' @aliases psdcore.default
#' @export
psdcore.default <- function(X.d, 
                            X.frq=NULL, 
                            ntaper=as.tapers(5), 
                            preproc=TRUE,
                            na.action = stats::na.fail,
                            plot=FALSE,
                            refresh=FALSE,
                            verbose=FALSE,
                            ndecimate,
                            ...
                           ) {
  
  if (!missing(ndecimate)) warning('Support for decimation has been removed.')
  
  # psd options
  ops <- getOption("psd.ops")
  stopifnot(!is.null(ops))
  evars <- ops[['names']]
  
  # named series
  series <- deparse(substitute(X.d))
  
  ## Convert to ts object
  if (!is.ts(X.d)){
    if (is.null(X.frq)){
      # make an assumption about the sampling rate
      X.frq <- 1
      if (verbose) message("\tsampling frequency not found -- taken to be\t", X.frq)
    }
    X.d <- if (X.frq > 0){
      # value represents sampling frequency
      stats::ts(X.d, frequency=X.frq)
    } else if (X.frq < 0){
      # value is sampling interval
      stats::ts(X.d, deltat=abs(X.frq))
    } else {
      stop("bad sampling information")
    }
  }
  # Smokey says: only you can stop NA fires
  X.d <- na.action(X.d)
  
  # Refresh sampling rate, and get Nyquist frequency, tapers, and status
  X.frq <- stats::frequency(X.d)
  Nyq <- X.frq / 2
  len_tapseq <- length(ntaper)
  single.taper.arg <- len_tapseq == 1
  
  # get a clean environment
  if ( single.taper.arg | refresh ) psd_envRefresh(verbose=verbose)
  
  #  only one variable in the env (init) means it hasn't been added to yet
  is.fresh <- length(psd_envStatus()[['listing']]) == 1
  
  ### Return complex fft, and initialize other things as necessary
  fftz <- if ( is.fresh ){
    #
    # initialize fft and other things, since this usually means a first run
    #
    # original series
    n.o <- psd_envAssignGet(evars[['n.orig']], length(X.d))
    X <- psd_envAssignGet(evars[['series.orig']], {
      if (preproc){
        # TODO: option for fast-detrend only, assign preproc flag in env (for plotting later)
        prewhiten(X.d, AR.max=0L, detrend=TRUE, plot=FALSE, verbose=verbose)$prew_lm
      } else {
        X.d
      }
    })
    
    # Force series to be even in length (modulo division)
    n.e <- psd_envAssignGet(evars[['n.even']], modulo_floor(n.o))
    X.even <- ts(X[seq_len(n.e)], frequency=X.frq)
    X.even <- psd_envAssignGet(evars[['series.even']], X.even)
    stopifnot(is.ts(X.even))
    
    # half length of even series
    nhalf <- psd_envAssignGet(evars[['n.even.half']], n.e/2)
    
    # variance of even series
    varx <- psd_envAssignGet(evars[['var.even']], drop(stats::var(X.even)))
    
    # create uniform tapers
    kseq <- psd_envAssignGet(evars[['last.taper']], {
      if (len_tapseq == 1){
        rep.int(ntaper, nhalf+1)
      } else {
        tmptap <- ntaper[seq_len(nhalf+1)] # if length < nhalf + 1 the remnants will be NA
        tmptap[is.na(tmptap)] <- pmin(ntaper[len_tapseq], ops[['tapmin']])
        tmptap
      }
    })
    
    ## zero pad and take double-length fft
    padded <- as.numeric(c(X.even, rep.int(0, n.e)))
    
    ## Calculate discrete Fourier tranform
    #   Note fftw is faster for very long series but we are
    #   using stats::fft until fftw is reliably built on CRAN
    padded.fft <- psd_envAssignGet(evars[['fft.padded']], stats::fft(padded))
    
    psd_envAssignGet(evars[['fft']], padded.fft)
  
  } else {
    
    if (verbose) warning("Working environment *not* refreshed. Results may be bogus.")
    
    X <- psd_envGet(evars[['series.even']])
    n.e <- psd_envGet(evars[['n.even']])
    nhalf <- psd_envGet(evars[['n.even.half']])
    varx <- psd_envGet(evars[['var.even']])
    kseq <- ntaper
    
    psd_envGet(evars[['fft']])
    
  }
  
  ### Switch: multitaper if TRUE, periodogram if FALSE
  do.mt <- if (single.taper.arg){
    ifelse(ntaper > 0, TRUE, FALSE)
  } else {
    ifelse(all(kseq > 0), TRUE, FALSE)
  }
  
  ###  Select frequencies for PSD evaluation
  f <- base::seq.int(0, nhalf, by=1)
  nfreq <- length(f)
  
  ###  Calculate the (un-normalized) PSD
  PSD <- psd_envAssignGet(evars[["last.psdcore"]], {
    
    if (do.mt){
      #
      # with sine tapers
      #
      if (verbose) message("\testimating multitaper psd")
      
      ## resample fft with taper sequence and quadratic weighting
      # ( this is where the majority of the computational work is )
      kseq <- as.integer(kseq) 
      reff <- resample_fft_rcpp(fftz, kseq, verbose=verbose)

      # return a valid resampled fft or stop
	    if (inherits(reff,'try-error')){
      	stop("Could not resample fft... inspect with psd_envGet(",evars[['fft']],"), etc.")
      } else {
        reff[['psd']]
      }
      
    } else {
      #
      # or with the simple cosine-tapered result
      #
      if (verbose) message("raw periodogram")
      
      Xfft <- psd_envGet(evars[['fft']])
      ff <- Xfft[seq_len(nfreq)]
      N0. <- psd_envGet(evars[['n.even']])
      
      # if the user wants it, then by all means
      # ... but force a warning on them
      warning("Careful interpreting raw-periodogram results!")
      base::abs(ff * base::Conj(ff)) / N0.
      
    }
  })
  
  # should not be complex at this point!
  stopifnot(!is.complex(PSD))
  
  # there should not be any bad values here!
  pNAs <- is.na(PSD)
  if (any(pNAs)) warning("NA psd estimates?!")
  
  ## Nyquist frequencies
  npsd <- length(PSD)
  frq <- as.numeric(base::seq.int(0, Nyq, length.out=npsd))
  
  ## Update tapers for consistency
  kseq <- as.tapers(if (do.mt){	
  	reff[['k.capped']]
  } else{ 
  	kseq
  })

  ## Normalize and convert to one-sided PSD
  #
  # ( using the trapezoidal rule, the principal being that the
  # integrated spectrum should be equal to the variance of the signal )
  #
  trap.area <- base::sum(PSD, na.rm=TRUE) - mean(PSD[c(1,npsd)], na.rm=TRUE)
  area.var.ratio <- varx / trap.area
  PSD <- 2 * PSD * area.var.ratio * nhalf
  
  ## Assemble final results
  mtap <- max(kseq, na.rm=TRUE)
  PSD.out <- list(freq = as.numeric(frq), 
                  spec = as.numeric(PSD), 
                  coh = NULL, 
                  phase = NULL, 
                  kernel = NULL, 
                  # must be a scalar for plot.spec to give conf ints:
                  df = 2 * mtap, # 2 DOF per taper, Percival and Walden eqn (370b)
                  numfreq = npsd,    
                  ## bandwidth
                  # http://biomet.oxfordjournals.org/content/82/1/201.full.pdf
                  # half-width W = (K + 1)/{2(N + 1)}
                  # effective bandwidth ~ 2 W (accurate for many spectral windows)
                  bandwidth = (mtap + 1) / nhalf, 
                  n.used = psd_envGet(evars[['n.even']]), 
                  orig.n = psd_envGet(evars[['n.orig']]), 
                  series = series, 
                  snames = colnames(X), 
                  method = sprintf("sine multitaper"), 
                  taper = kseq, 
                  pad = TRUE, # always!
                  detrend = preproc, # always true?
                  demean = preproc,
                  timebp = as.numeric(kseq/2),
                  nyquist.frequency = Nyq
  )
  class(PSD.out) <- c("amt","spec")
  if (plot) pgram_compare(PSD.out, ...)
  return(invisible(PSD.out))
}

#' Compare multitaper spectrum with cosine-tapered periodogram
#' 
#' @description
#' Plot the results of \code{\link{psdcore}} against the results of 
#' \code{\link[stats]{spec.pgram}}
#' 
#' @name pgram_compare
#' @export
#' 
#' @param x a single \code{\link{psdcore}} object
#' @param f numeric; the frequency range to plot; optional: if not given the program will show the entire band.
#' @param X object used to create \code{x}; optional: if not given the program will
#' try and access the last copy in the environment. An attempt is made to coerce to an object of class \code{'ts'}.
#' @param log.freq logical; should frequencies be transformed with \code{\link{log10}}?  
#' Note that if \code{f} is given, the values should not already be transformed.
#' @param db.spec logical; should the spectrum estimates be converted to decibels with \code{\link{dB}}?
#' @param taper numeric; specifies the proportion of data to taper for the cosine periodogram. 
#' @param ... additional parameters (currently unused)
#' 
#' @return A list with the cosine-tapered estimates and the adaptive estimates, invisibly.
#' @examples
#' set.seed(1234)
#' X <- rnorm(1e3)
#' 
#' # multitaper spectrum
#' p <- psdcore(X, ntaper=10)
#' 
#' # how does it compare to a single-cosine tapered spectrum?
#' pgram_compare(p)
#' 
#' # or in a certain band
#' pgram_compare(p, c(0.1,0.4))
#' 
#' # linear frequencies
#' pgram_compare(p, c(0.1,0.4), log.freq = FALSE)
pgram_compare <- function(x, ...) UseMethod("pgram_compare")

#' @rdname pgram_compare
#' @aliases pgram_compare.amt
#' @export
pgram_compare.amt <- function(x, f=NULL, X=NULL, log.freq=TRUE, db.spec=TRUE, taper=0.2, ...){
  
  P. <- x
  f. <- P.[['freq']]
  s. <- P.[['spec']]
  t. <- P.[['taper']]
  
  flims <- if (!is.null(f)){
    if (log.freq) {
      log10(f)
    } else {
      f
    }
  } else {
    if (log.freq) {
      log10(range(f.[f.>0]))
    } else {
      range(f.)
    }
  }
  if (any(is.infinite(flims))) flims <- NULL
  
  detrend <- P.[['detrend']]
  demean <- P.[['demean']]
  
  ops <- getOption("psd.ops")
  stopifnot(!is.null(ops))
  evars <- ops[['names']]
  
  if (is.null(X)) X <- psd_envGet(evars[['series.even']])
  stopifnot(!is.null(X))
  X <- as.ts(X)
  fsamp <- frequency(X)
  
  # Calculate (single) cosine taper psd
  tc. <- taper
  Pc. <- spec.pgram(X, log="no", pad=1, taper=tc., detrend=detrend, demean=detrend, plot=FALSE)
  # frequencies are appropriate,
  # but spectrum is normed for double-sided whereas psd single-sided; hence,
  # factor of 2
  Pc. <- normalize(Pc., fsamp, "spectrum", verbose=FALSE)
  fc. <- Pc.[['freq']]
  sc. <- Pc.[['spec']]
  
  flab <- if (log.freq){
    f. <- log10(f.)
    fc. <- log10(fc.)
    expression(log[10] ~~ "frequency")
  } else {
    expression('frequency')
  }
  
  slab <- if (db.spec){
    s. <- dB(s./fsamp) # normalization adjustment
    sc. <- dB(sc.)
    expression("dB rel." ~~ "units"**2 %*% delta * t)
  } else {
    expression(delta ~ "units"**2 %*% delta * t)
  }
  slims <- range(pretty(c(s.,sc.)))
  
  ## Start plotting
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  
  #layout(matrix(c(1,2), ncol=1), c(1,2))
  layout(matrix(c(1,1,2,2,3,4), 3, 2, byrow = TRUE))
  par(mar=c(2, 3.5, 2.3, 1.2), oma=c(2,2,2,1))
  par(las=1, tcl = -.3, mgp=c(2.2, 0.4, 0))
  par(cex=0.8)
  
  ## Periodogram result first
  plot(fc., sc., 
       col="red", type="l",
       xaxs="i", xlim=flims, 
       yaxs="i", ylim=slims, 
       xlab="", 
       ylab=slab,
       main="")
  # and adaptive result overlain
  lines(f., s., type="l")
  legend("bottomleft",
         c(sprintf("spec.pgram (%i%% cosine taper)", 100*tc.), 
           sprintf("psdcore (%i - %i tapers)", min(t.), max(t.)) ), 
         col=c("red","black"), lty=1, lwd=2, cex=0.9)
  mtext(flab, side=1, line=2)
  mtext("Spectral density estimates", adj=0, line=0.2, font=4)
  
  ## Tapers
  if (is.tapers(t.)){
    plot(t., f., xlim=flims, xaxs="i")
  } else {
    plot(f., t., type="h", xlim=flims, xaxs="i")
  }
  mtext("Sine-tapers", adj=0, line=0.2, font=4)
  #mtext(flab, side=1, line=2)
    
  ## Original series
  main.pre <- ifelse(detrend | demean, "Modified e", "E")
  main.post <- sprintf("(dem. %s detr. %s)", detrend, demean) 
  plot(X, type="l", ylab="units", xlab="", xaxs="i", main="")
  mtext(paste0(main.pre, "ven-length series"), line=1, font=4)
  mtext(main.post, cex=0.6)
  mtext(expression("time"), side=1, line=1.7)
  
  ## autocorrelation
  acf(X, main="")
  mtext(expression("lag, time"), side=1, line=1.7)
  mtext("Auto-correlation function", line=1, font=4)
  
  title("Sine-multitaper PSD vs. Tapered Periodogram", outer=TRUE, line=0, cex.main=1.5)
  return(invisible(list(spec=Pc., amt=P.)))
  
}
