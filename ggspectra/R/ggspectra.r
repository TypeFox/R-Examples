#' @details
#' Pakage \code{ggspectra} is a package with extensions to ggplot2 for plotting
#' radiation spectra. It is of little use for other types of data as most
#' functions depend on x aesthetic being mapped to a variable containing
#' wavelength values expressed in nanometres. It works well together with
#' some other extensions to package 'ggplot2' such as pacakegs 'ggrepel' and
#' 'cowplot'. This package is very tightly dependent on package
#' \code{\link[photobiology]{photobiology}}.
#'
#' @references
#' \code{ggplot2} web site at \url{http://ggplot2.org/}\cr
#' \code{ggplot2} source code at \url{https://github.com/hadley/ggplot2}\cr
#' Function \code{multiplot} from \url{http://www.cookbook-r.com/}
#'
#' @author Pedro J. Aphalo
#'
#' @import photobiology photobiologyWavebands ggplot2
#' @importFrom graphics plot
#'
#' @note
#' This package makes use of the new features of 'ggplot2' 2.0.0 that make
#' writing this kind of extensions really easy and is consequently not
#' compatible with earlier versions of 'ggplot2'.
#'
#' @examples
#' library(ggplot2)
#' library(photobiology)
#' library(photobiologyWavebands)
#'
#' # maximum
#' ggplot(sun.spct, aes(w.length, s.e.irrad)) + geom_line() +
#' stat_peaks(span = NULL)
#'
#' ggplot(sun.spct, aes(w.length, s.e.irrad)) + geom_line() +
#' stat_peaks(span = 21, geom = "text")
#'
#' ggplot(sun.spct, aes(w.length, s.e.irrad)) + geom_line() +
#'   stat_valleys(span = 21, geom = "text")
#'
#' ggplot(sun.spct, aes(w.length, s.e.irrad)) + geom_line() +
#'   stat_peaks(span = 21, geom = "point", colour = "red") +
#'   stat_peaks(span = 51, geom = "text", colour = "red", vjust = -0.3,
#'              label.fmt = "%3.0f nm")
#'
#' ggplot(sun.spct, aes(w.length, s.e.irrad)) + geom_line() +
#'   stat_color() + scale_color_identity()
#'
#' plot(sun.spct)
#' plot(polyester.spct, UV_bands(), range = UV())
#'
"_PACKAGE"
