#' CDOM absorption data.
#'
#' Simple absorption spectra used to test package's functions.
#'
#' \itemize{
#'   \item wavelength.  Wavelengths used for measurements (190-900 nm.)
#'   \item Absorption
#' }
#'
#' @import ggplot2
#' @import tidyr
#' @docType data
#' @keywords datasets
#' @name spectra
#' @usage data(spectra)
#' @format A data frame with 711 rows and 26 variables
#' @examples
#' library(ggplot2)
#' library(tidyr)
#' data("spectra")
#'
#' spectra <- gather(spectra, sample, absorption, -wavelength)
#'
#' ggplot(spectra, aes(x = wavelength, y = absorption, group = sample)) +
#'  geom_line(size = 0.1)
NULL
