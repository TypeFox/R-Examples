#' Calculate a normalized index.
#'
#' This function returns a normalized difference index value for an arbitrary
#' pair of wavebands.
#'
#' @param spct an R object
#' @param plus.w.band waveband objects The waveband determine the
#'   region of the spectrum used in the calculations
#' @param minus.w.band waveband objects The waveband determine the
#'   region of the spectrum used in the calculations
#' @param f function used for integration taking spct as first argument and a
#'   list of wavebands as second argument.
#' @param ... additional arguments passed to f
#'
#' @return A numeric value for the index
#'
#' @export
#'
#' @note \code{f} is most frequently \code{\link{reflectance}}, but also
#'   \code{\link{transmittance}}, or even \code{\link{absorbance}},
#'   \code{\link{response}}, \code{\link{irradiance}} or a user-defined function
#'   can be used if there is a good reason for it. In every case \code{spct}
#'   should be of the class expected by \code{f}. When using two wavebands of
#'   different widths do consider passing to \code{f} a suitable \code{quantity}
#'   argument. Wavebands can describe weighting functions if desired.
#'
#' @export
#'
normalized_diff_ind <- function(spct, plus.w.band, minus.w.band, f, ...) {
  x <- as.numeric(f(spct, list(plus.w.band, minus.w.band), ...))
  z <- (x[1] - x[2]) / (x[1] + x[2])
  name <- paste("NDI ", as.character(substitute(f)), " [",
                sub("range.", "", labels(plus.w.band)$label), "] - [",
                sub("range.", "", labels(minus.w.band)$label), "]")
  names(z) <- name
  z
}

