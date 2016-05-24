#' Spectrum to rgb color conversion
#'
#' Calculates rgb values from spectra based on human color matching functions
#' (CMF) or chromaticity coordinates (CC). A CMF takes into account luminous
#' sensitivity, while a CC only the color hue. This function, in contrast to
#' that in package pavo does not normalize the values to equal luminosity, so
#' using a CMF as input gives the expected result. Another difference is that it
#' allows the user to choose the chromaticity data to be used. The data used by
#' default is different, and it corresponds to the whole range of CIE standard,
#' rather than the reduced range 400 nm to 700 nm. The wavelength limits are not
#' hard coded, so the function could be used to simulate vision in other
#' organisms as long as pseudo CMF or CC data are available for the simulation.
#'
#' @param w.length numeric array of wavelengths (nm)
#' @param s.e.irrad numeric array of spectral irradiance values
#' @param sens a chroma_spct object with variables w.length, x, y, and z, giving
#'   the CC or CMF definition (default is the proposed human CMF according to
#'   CIE 2006.)
#' @param color.name character string for naming the rgb color definition
#' @param check logical indicating whether to check or not spectral data
#'
#' @return A color defined using \code{\link[grDevices]{rgb}}. The numeric
#'   values of the RGB components can be obtained using function
#'   \code{\link[grDevices]{col2rgb}}.
#'
#' @export
#'
#' @examples
#' my.color <-
#'     with(sun.data, s_e_irrad2rgb(w.length, s.e.irrad, color.name="sunWhite"))
#' col2rgb(my.color)
#'
#' @note Very heavily modified from Chad Eliason's
#'   \email{cme16@@zips.uakron.edu} spec2rgb function in package \code{Pavo}.
#' @references CIE(1932). Commission Internationale de l'Eclairage Proceedings,
#'   1931. Cambridge: Cambridge University Press.
#' @references Color matching functions obtained from Colour and Vision Research
#' Laboratory online data respository at \url{http://www.cvrl.org/}.
#'
#' @seealso \url{http://www.cs.rit.edu/~ncs/color/t_spectr.html}.
#'
#' @family color functions
#'
s_e_irrad2rgb <- function(w.length, s.e.irrad,
                          sens=photobiology::ciexyzCMF2.spct,
                          color.name=NULL, check=TRUE) {
  low.limit <- min(sens$w.length)
  high.limit <- max(sens$w.length)
  if (single_wl <- length(w.length) == 1) {
    if (w.length < low.limit || w.length > high.limit) {
      return(grDevices::rgb(0, 0, 0, names=color.name))
    } else {
      s.e.irrad = 1.0
    }
  } else {
    if (check && !check_spectrum(w.length, s.e.irrad)) {
      return(NA)
    }
}
# if we have a spectrum we will expand and fill with zeros when needed

if (!single_wl) {
  if ((max(w.length) <= low.limit) || (min(w.length) >= high.limit)) {
    return("black")
  }
  sens$s.e.irrad <- interpolate_spectrum(w.length, s.e.irrad, sens$w.length, fill=0.0)
  sens$s.e.irrad.norm <- with(sens, s.e.irrad / integrate_xy(w.length, s.e.irrad))

  X <- with(sens, integrate_xy(w.length, s.e.irrad.norm * x))
  Y <- with(sens, integrate_xy(w.length, s.e.irrad.norm * y))
  Z <- with(sens, integrate_xy(w.length, s.e.irrad.norm * z))
} else {
  X <- stats::approx(sens$w.length, sens$x, w.length)$y
  Y <- stats::approx(sens$w.length, sens$y, w.length)$y
  Z <- stats::approx(sens$w.length, sens$z, w.length)$y
}

XYZ <- rbind(X, Y, Z)

xyzmat <- rbind(c(3.240479, -1.537150, -0.498535),
                c(-0.969256, 1.875992, 0.041556),
                c(0.055648, -0.204043, 1.057311))

rgb1 <- xyzmat %*% as.matrix(XYZ)

# print(rgb1)

# not all colours can be represented in the RGB space, so unrepresentable
# colours are converted by brute force into representable colours
rgb1[rgb1 < 0] <- 0
rgb1[rgb1 > 1] <- 1

if (anyNA(rgb1[ , 1])) {
  warning("NA in rgb values, returning 'black'")
  return("black")
}

rgb.color <- grDevices::rgb(red=rgb1[1,1], green=rgb1[2,1], blue=rgb1[3,1], names=color.name)

rgb.color

}
