#' Wavelength to rgb color conversion
#'
#' Calculates rgb values from spectra based on human color matching functions
#'
#' @param w.length numeric Vector of wavelengths (nm)
#' @param sens chroma_spct Used as chromaticity definition
#' @param color.name character Used for naming the rgb color definition
#'
#' @return A vector of colors defined using \code{rgb()}. The numeric values of
#'   the RGB components can be obtained using function \code{col2rgb()}.
#'
#' @export
#' @examples
#' col2rgb(w_length2rgb(580))
#' col2rgb(w_length2rgb(c(400, 500, 600, 700)))
#' col2rgb(w_length2rgb(c(400, 500, 600, 700), color.name=c("a","b","c","d")))
#' col2rgb(w_length2rgb(c(400, 500, 600, 700), color.name="a"))
#'
#' @family color functions
#'
w_length2rgb <- function(w.length,
                         sens = photobiology::ciexyzCMF2.spct,
                         color.name = NULL) {
  len.wl <- length(w.length)
  generate.names <- is.null(color.name)
  if (!generate.names) {
    len.col <- length(color.name)
    if (len.col == 1L) {
      color.names <- rep(color.name[1], length.out = len.wl)
    } else if (len.col < len.wl) {
      warning("color.name argument shorter than w.length argument.")
      color.names <- color.name
    } else {
      color.names <- color.name
    }
  } else {
    color.names <-NULL
  }
   colors <- NULL
   for (i in 1:len.wl) {
     colors[i] <-  s_e_irrad2rgb(w.length[i], 1.0, sens=sens, check = FALSE)
   }
  if (generate.names) {
    color.names <- paste("wl", as.character(round(w.length, 1)), "nm", sep=".")
  }
  names(colors) <- color.names
  return(colors)
}
