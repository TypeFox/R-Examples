#' Wavelength range to rgb color conversion
#'
#' Calculates rgb values from spectra based on human color matching functions
#'
#' @param w.length numeric Vector of wavelengths (nm) of length 2. If longer,
#'   its range is used.
#' @param sens chroma_spct Used as the chromaticity definition
#' @param color.name character Used for naming the rgb color definition
#'
#' @return A vector of colors defined using \code{rgb()}. The numeric values of
#'   the RGB components can be obtained using function \code{col2rgb()}.
#'
#' @export
#' @examples
#' col2rgb(w_length_range2rgb(c(500,600)))
#' col2rgb(w_length_range2rgb(550))
#' col2rgb(w_length_range2rgb(500:600))
#'
#' @family color functions
#'
w_length_range2rgb <- function(w.length,
                               sens=photobiology::ciexyzCMF2.spct,
                               color.name=NULL) {
  if (is.null(w.length) || !is.numeric(w.length)) {
    warning("Bad wlength input, must be numeric")
    return("black")
  }
  w.length <- unique(sort(w.length))
  len <- length(w.length)
  if (len == 1 ) {
    message("Calculating RGB values for monochromatic light.")
    return(w_length2rgb(w.length, sens, color.name))
  } else if (len > 2) {
    message("Using only extreme wavelength values.")
    w.length <- range(w.length)
  } else if (len < 1) {
    stop("Bad assertion!")
  }
  num.values <- min(5L, ceiling(spread(w.length)))
  w.length.values <- seq(w.length[1], w.length[2], length.out = num.values)
  s.e.irrad.values <- rep(1.0, length.out = num.values)
  color <-  s_e_irrad2rgb(w.length.values, s.e.irrad.values, sens=sens,
                                color.name=ifelse(is.null(color.name),
                              paste(as.character(w.length[1]), "-", as.character(w.length[2]), " nm", sep=""),
                           color.name), check=FALSE)
  return(color)
}
