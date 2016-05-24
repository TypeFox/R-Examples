#' @details Package 'photobiology' is at the core of a suite of packages for
#'   analysis and plotting of data relevant to photobiology (described at
#'   \url{http://www.r4photobiology.info/}). The accompanying packages (under
#'   development) provide data and definitions that are to a large extent
#'   application-area specific while the functions in the present package are
#'   widely useful in photobiology and radiation quantification in geophysics
#'   and meteorology. Package 'photobiology' has its main focus in the
#'   characterization of the light environment in a biologically relevant manner
#'   and in the manipulation of spectral data to simulate photo-physical,
#'   photo-chemical and photo-biological interactions and reponses. The focus of
#'   package 'pavo' (Maia et al., 2003) is in colour perception by animals and
#'   assessment of animal coloration. In spite of the different focus, there is
#'   some degree of overlap.
#'
#' @note Code for some of the astronomical calculations has been adapted from
#'   that in package 'pavo'.
#'
#' @references
#' Aphalo, P. J., Albert, A., Björn, L. O., McLeod, A. R., Robson, T. M.,
#' Rosenqvist, E. (Eds.). (2012). Beyond the Visible: A handbook of best
#' practice in plant UV photobiology (1st ed., p. xxx + 174). Helsinki:
#' University of Helsinki, Department of Biosciences, Division of Plant Biology.
#' ISBN 978-952-10-8363-1 (PDF), 978-952-10-8362-4 (paperback). Open access PDF
#' download available at http://hdl.handle.net/10138/37558
#'
#' Maia, R., Eliason, C. M., Bitton, P. P., Doucet, S. M., Shawkey, M. D. (2013)
#' pavo: an R package for the analysis, visualization and organization of
#' spectral data. Methods in Ecology and Evolution, 4(10):906-913. doi:
#' 10.1111/2041-210X.12069
#'
#' @section Acknowledgements: This work was funded by the Academy of Finland
#'   (decision 252548). COST Action FA9604 'UV4Growth' facilitated discussions
#'   and exchanges of ideas that lead to the development of this package. The
#'   contributins of Andy McLeod, Lars Olof Björn, Nigel Paul, Lasse Ylianttila,
#'   T. Matthew Robson and Titta Kotilainen were specially significant.
#'   Tutorials by Hadley Wickham comments on my presentation at UseR!2015
#'   allowed me to significantly improving the coding and functionality.
#'
#' @examples
#' # irradiance of the whole spectrum
#' irrad(sun.spct)
#' # photon irradiance 400 nm to 700 nm
#' q_irrad(sun.spct, waveband(c(400,700)))
#' # energy irradiance 400 nm to 700 nm
#' e_irrad(sun.spct, waveband(c(400,700)))
#' # simulating the effect of a filter on solar irradiance
#' e_irrad(sun.spct * yellow_gel.spct, waveband(c(400,500)))
#' e_irrad(sun.spct * yellow_gel.spct, waveband(c(500,700)))
#' # daylength
#' sunrise_time(lubridate::today(tzone = "EET"), tz = "EET",
#'              lat = 60, lon = 25, unit.out = "hour")
#' day_length(lubridate::today(tzone = "EET"), tz = "EET",
#'               lat = 60, lon = 25, unit.out = "hour")
#' # colour as seen by humans
#' color(sun.spct)
#' color(sun.spct * yellow_gel.spct)
#' # filter transmittance
#' transmittance(yellow_gel.spct)
#' transmittance(yellow_gel.spct, waveband(c(400,500)))
#' transmittance(yellow_gel.spct, waveband(c(500,700)))
"_PACKAGE"
