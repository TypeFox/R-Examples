#' CDOM absorbance data.
#'
#' Simple absorbance spectra used to test package's functions.
#'
#' \itemize{
#'   \item wavelength.  Wavelengths used for measurements (190-900 nm.)
#'   \item absorbance
#' }
#' @docType data
#' @keywords datasets
#' @name absorbance
#' @usage data("absorbance")
#' @format A data frame with 711 rows and 4 variables
#' @examples
#'
#' data("absorbance")
#'
#' plot(absorbance$wavelength, absorbance$sample1, type = "l",
#' xlab = "Wavelengths", ylab = "Absorbance per meter")
#' lines(absorbance$wavelength, absorbance$sample2, col = "blue")
#' lines(absorbance$wavelength, absorbance$sample3, col = "red")

NULL
