#' @title Estimates spectrum of multivariate time series
#' @name mvspectrum
#' @importFrom sapa SDF
#' @importFrom ifultools checkScalarType
#' @description
#' The spectrum of a multivariate time series is a matrix-valued function of the 
#' frequency \eqn{\lambda \in [-\pi, \pi]}, which is symmetric/Hermitian around 
#' \eqn{\lambda = 0}.
#' 
#' \code{mvspectrum} estimates it and returns a 3D array of dimension 
#' \eqn{num.freqs \times K \times K}.  Since the spectrum is symmetric/Hermitian around
#'  \eqn{\lambda = 0} it is sufficient to store only positive frequencies.  
#' In the implementation in this package we thus usually 
#' consider only positive frequencies (omitting \eqn{0}); \code{num.freqs} refers
#' to the number of positive frequencies only.
#' 
#' @inheritParams common-arguments
#' @param method string; method for spectrum estimation; see \code{method} argument in
#' \code{\link[sapa]{SDF}} (in the \pkg{sapa} package); use 
#' \code{"mvspec"} to use \code{\link[astsa]{mvspec}} (\pkg{astsa} package); or
#' use \code{"pgram"} to use \code{\link[stats]{spec.pgram}}.
#' @param normalize logical; if \code{TRUE} the spectrum will be normalized (see 
#' Value below for details).
#' @param \dots additional arguments passed to \code{\link[sapa]{SDF}} or 
#' \code{\link[astsa]{mvspec}} (e.g., \code{taper})
#' @return 
#' \code{mvspectrum} returns a 3D array of dimension \eqn{num.freqs \times K \times K}, where
#' \itemize{
#'  \item num.freqs is the number of frequencies
#'  \item K is the number of series (columns in \code{series}).
#' }
#' Note that it also has an attribute \code{"normalized"} which is
#' \code{FALSE} if \code{normalize = FALSE}; otherwise \code{TRUE}.
#' See \code{normalize_mvspectrum} for details.
#' @references 
#' See References in \code{\link[stats]{spectrum}}, \code{\link[sapa]{SDF}}, 
#' \code{\link[astsa]{mvspec}}.
#' @keywords ts
#' @export
#' @examples
#' 
#' set.seed(1)
#' XX <- cbind(rnorm(100), arima.sim(n = 100, list(ar = 0.9)))
#' ss3d <- mvspectrum(XX)
#' dim(ss3d)
#' 
#' ss3d[2,,] # at omega_1; in general complex-valued, but Hermitian
#' identical(ss3d[2,,], Conj(t(ss3d[2,,]))) # is Hermitian
#' 
#' ss <- mvspectrum(XX[, 1], smoothing = TRUE)
#' 
#' \dontrun{
#'   mvspectrum(XX, normalize = TRUE)
#' }
#' ss <- mvspectrum(whiten(XX)$U, normalize = TRUE)
#' 
#' xx <- scale(rnorm(100), center = TRUE, scale = FALSE)
#' var(xx)
#' sum(mvspectrum(xx, normalize = FALSE, method = "direct")) * 2
#' sum(mvspectrum(xx, normalize = FALSE, method = "wosa")) * 2
#' 
#' 

mvspectrum <- function(series, 
                       method = 
                         c("pgram", "multitaper", "direct", "wosa", 
                           "mvspec", "ar"),
                       normalize = FALSE, smoothing = FALSE, ...) {
  
  method <- match.arg(method)
  
  addl.args <- list(...)
  if ("spectrum.control" %in% addl.args) {
    stop("'spectrum.control' is not a valid argument of 'mvspectrum'.",
         "Use method = spectrum.control$method and ",
         "smoothing = spectrum.control$smoothing directly instead.")
  } 
  if ("entropy.control" %in% addl.args) {
    stop("'entropy.control' is not a valid argument of 'mvspectrum'.")
  }
  
  series <- as.matrix(series)
  num.series <- ncol(series)
  num.obs <- nrow(series)
  
  if (num.series > 1 && smoothing) {
    warning("Smoothing can only be set to TRUE for univariate time series.\n",
            "Will be ignored for multivariate time series.")
  }
  
  if (normalize) {
    # if the data is not uncorrelated with unit variance, 
    # then dont allow to normalize
    series <- check_whitened(series)
  }
  
  if (method == "mvspec") {
    stopifnot(requireNamespace("astsa", quietly = TRUE))
    out <- .mvspec2mvspectrum(astsa::mvspec(series, plot = FALSE, 
                                            detrend = FALSE, fast = FALSE, 
                                            ...))
  } else if (method == "ar") {
    stopifnot(num.series == 1)
    out <- spec.ar(c(series), method = "burg", plot = FALSE,
                   n.freq = ceiling(length(series) / 2 + 1))$spec[-1]
  } else if (method %in% c("wosa", "multitaper", "direct")) {
    stopifnot(requireNamespace("sapa", quietly = TRUE)) 
    out <- .SDF2mvspectrum(sdf.output = sapa::SDF(series, method = method, 
                                                  recenter = TRUE, 
                                                  npad = num.obs, ...))
    if (num.series == 1) {
      out <- array(out, dim = c(length(out), 1, 1))
    } else {
      out <- as.array(out)
    }
  } else if (method == "pgram") {
    out <- mvpgram(series)
  } else {
    stop("Method '", method, "' is not implemented.\n", 
         "Please specify other method.")
  }
  
  if (smoothing && num.series == 1) {
    if (requireNamespace("mgcv", quietly = TRUE)) {
      mod.tmp <- mgcv::gam(c(out) ~ s(seq_along(out)),
                           family = Gamma(link = "log"), method = "REML")
      dim.out <- dim(out)
      # dispersion = 1 for an exponential distribution
      out <- predict(mod.tmp, type = "response", dispersion = 1)
      dim(out) <- dim.out
    } else {
      warning("Spectrum was not smoothed since this requires the 'mgcv' package.\n", 
              "Please install it.")
    }
  }
  
  # divide by number of observations  = 2 * number of positive frequencies
  # This takes case of 
  #   - divide by 2 since it is symmetric around 0
  #   - divide by number of frequencies so that the sum over spectrum = variance
  out <- out / num.obs
 
  # add frequencies from (0, pi]; remove '0' frequency
  if (num.series == 1) {
    num.freqs.in.0.pi <- length(out)
  } else {
    num.freqs.in.0.pi <- nrow(out)
  }
  attr(out, "frequency") <- 
    seq(1, num.freqs.in.0.pi, by = 1) / (2 * num.freqs.in.0.pi + 1) * (2 * pi)
  attr(out, "normalized") <- FALSE

  if (normalize) {
    out <- normalize_mvspectrum(out)
  }
  
  class(out) <- "mvspectrum"
  invisible(out)
} 
#' @rdname mvspectrum
#' @description
#' \code{normalize_mvspectrum} normalizes the spectrum such that
#' it adds up to \eqn{0.5} over all positive frequencies (by symmetry it will 
#' add up to 1 over the whole range -- thus the name \emph{normalize}). 
#' 
#' For a \eqn{K}-dimensional time series it adds
#' up to a Hermitian \eqn{K \times K} matrix with 0.5 in the diagonal and
#' imaginary elements (real parts equal to \eqn{0}) in the off-diagonal. 
#' Since it is Hermitian the mvspectrum will add up to the identity matrix
#' over the whole range of frequencies, since the off-diagonal elements
#' are purely imaginary (real part equals 0) and thus add up to 0.
#' 
#' @details
#' For an orthonormal time series \eqn{\mathbf{U}_t} the raw periodogram adds up 
#' to \eqn{I_K} 
#' over all (negative and positive) frequencies.  Since we only consider
#' positive frequencies, the normalized multivariate spectrum should add up to
#' \eqn{0.5 \cdot I_K} plus a Hermitian imaginary matrix (which will add up to zero
#' when combined with its symmetric counterpart.)
#' As we often use non-parametric smoothing for less variance, the spectrum estimates
#' do not satisfy this identity exactly. \code{normalize_mvspectrum} thus adjust the 
#' estimates so they satisfy it again exactly. 
#'  
#' @keywords manip
#' @inheritParams common-arguments
#' @export
#' @return
#' \code{normalize_mvspectrum} returns a normalized spectrum over 
#' positive frequencies, which:
#' \describe{
#'   \item{univariate:}{adds up to \eqn{0.5},}
#'   \item{multivariate:}{adds up to Hermitian \eqn{K \times K} matrix
#'   with 0.5 in the diagonal and purely imaginary elements in the off-diagonal.} 
#' }
#' 
#' @examples
#' xx <- scale(rnorm(100), center = TRUE, scale = FALSE)
#' ss <- mvspectrum(xx)
#' ss.n <- normalize_mvspectrum(ss)
#' sum(ss.n)
#' # multivariate
#' UU <- whiten(matrix(rnorm(40), ncol = 2))$U
#' S.U <- mvspectrum(UU, method = "wosa")
#' mvspectrum2wcov(normalize_mvspectrum(S.U))
#' 

normalize_mvspectrum <- function(mvspectrum.output) {
  
  if (is.null(dim(mvspectrum.output))) {
    num.series <- 1
    num.freqs <- length(mvspectrum.output)
  } else {
    num.series <- ncol(mvspectrum.output) 
    num.freqs <- nrow(mvspectrum.output)
  }
  
  if (num.series == 1) {
    # just one dimensional series
    f3D <- mvspectrum.output / sum(mvspectrum.output)
    num.series <- 1
  } else {
    # note that this multiplies by 2 already; as we later us this
    # ad a divisor (and thus divide by 2), we have to multiply by 2 again
    # below
    cov.est.spectrum <- mvspectrum2wcov(mvspectrum.output)
    ss.sqrt.inv <- sqrt_matrix(cov.est.spectrum, return.sqrt.only = FALSE,
                               symmetric = TRUE)$sqrt.inverse
    # here multiply by 2 again as we use the full left and right spectrum
    # to estimate the covariance matrix
    f3D <- array(t(apply(mvspectrum.output, 1, 
                         function(x) 2 * t(Conj(ss.sqrt.inv)) %*% x %*%
                           ss.sqrt.inv)),
                 dim(mvspectrum.output))
  }
  # divide by 2 since we only consider positive frequencies here 
  # (and it is symmetric around 0)
  attr(f3D, "normalized") <- TRUE
  attr(f3D, "frequency") <- attr(mvspectrum.output, "frequency")
  return(f3D / 2)
} 


#' @rdname mvspectrum
#' @export
#' @description
#' \code{check_mvspectrum_normalized} checks if the spectrum is normalized 
#' (see \code{\link{normalize_mvspectrum}} for the requirements).
#' @inheritParams common-arguments
#' @param check.attribute.only logical; if \code{TRUE} it checks the 
#' attribute only.  This is much faster (it just needs to look up one attribute
#' value), but it might not surface silent bugs.  For sake of performance
#' the package uses the attribute version by default.  However, for 
#' testing/debugging the full computational version can be used.
#' @return
#' \code{check_mvspectrum_normalized} throws an error if spectrum is not
#' normalized correctly.
#' 
check_mvspectrum_normalized <- function(f.U, check.attribute.only = TRUE) {
  
  if (check.attribute.only) {
    stopifnot(attr(f.U, "normalized") == TRUE)
  } else {
    if (is.null(dim(f.U))) {
      num.series <- 1
      sum.of.spectrum <- as.matrix(sum(f.U))
    } else {
      num.series <- dim(f.U)[2]
      sum.of.spectrum <- apply(f.U, 2:3, sum)
    }
  
    off.diag <- sum.of.spectrum
    diag(off.diag) <- 0
    
    equals <- list()
    equals[["add-to-0.5"]] <- all.equal(target = diag(0.5, num.series),
                                        current = Re(sum.of.spectrum),
                                        check.names = FALSE,
                                        check.attributes = FALSE)
    equals[["off-diagonal-zero"]] <- all.equal(target = matrix(0, num.series, num.series),
                                               current = Re(off.diag),
                                               check.names = FALSE,
                                               check.attributes = FALSE)
    equals[["Hermitian"]] <- all.equal(target = off.diag,
                                       current = Conj(t(off.diag)),
                                       check.names = FALSE,
                                       check.attributes = FALSE)
    
    errors <- c("add-to-0.5" = paste("A _normalized_ spectrum must add up to",
                                          "0.5 times identity matrix."),
                "off-diagonal-zero" = paste("Off diagonals must have zero reals."),
                "Hermitian" = paste("The normalized spectrum must be Hermitian."))
    
    ind.errors <- !sapply(equals, isTRUE)
    
    if (any(ind.errors) > 0) {
      for (ii in seq_along(ind.errors)) {
        if (ind.errors[ii]) {
          cat(names(equals)[ii], ": ", equals[[ii]], "\n", sep = "")
        }
      }
      stop("Spectral density must be normalized:\n \t ", 
           paste(errors[ind.errors], collapse = "\n \t ")) 
    }
  }
}


#' @rdname mvspectrum
#' @export
#' @description
#' \code{mvpgram} computes the multivariate periodogram estimate using
#' bare-bone multivariate fft (\code{\link[stats]{mvfft}}). Please use
#' \code{mvspectrum(..., method = 'pgram')} instead of \code{mvpgram} directly.
#' 
#' This function is merely included to have one method that does not
#' require the \pkg{astsa} nor the \pkg{sapa} R packages.  However, 
#' it is strongly encouraged to install either one of them to get (much)
#' better estimates.  See Details.
#' 
#' @details
#' \code{mvpgram} has no options for improving spectrum estimation whatsoever.
#' It thus yields very noisy (in fact, inconsistent) estimates of the 
#' multivariate spectrum \eqn{f_{\mathbf{X}}(\lambda)}. 
#' If you want to obtain better estimates then please use other \code{method}s in
#' \code{\link{mvspectrum}} (this is highly recommended to obtain more 
#' reasonable/stable estimates).

mvpgram <- function(series) {
  if (is.null(dim(series))) {
    num.obs <- length(series)
    num.series <- 1
  } else {
    num.obs <- nrow(series)
    num.series <- ncol(series)
  }
  
  num.freqs <- floor(num.obs / 2)
  if (num.series > 1) {
    series.fft <- mvfft(series)
  } else {
    series.fft <- as.matrix(fft(series), ncol = 1)
  }
  # Periodogram as in spec.pgram
  pgram.x <- array(NA, dim = c(num.obs, num.series, num.series))
  
  for (ii in seq_len(num.series)) {
    for (jj in seq_len(num.series)) {
      pgram.x[, ii, jj] <- series.fft[, ii] * Conj(series.fft[, jj]) / num.obs
    }
  }
  # remove 0 frequency and only take positive frequencies
  pgram.x <- pgram.x[1 + seq_len(num.freqs),,]
  if (num.series == 1) {
    pgram.x <- Re(pgram.x)
  }
  return(pgram.x)
}