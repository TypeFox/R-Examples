#' Setlow's action spectrum for DNA damage
#'
#' A dataset containing the wavelengths at a 0.1 nm interval. Tabulated values
#' for Setlow's naked DNA damage action spectrum as used in the TUV model.
#'
#' The variables are as follows:
#'
#' \itemize{
#'   \item w.length (nm)
#'   \item s.e.response
#' }
#'
#' @docType data
#' @keywords datasets
#' @format A response.spct object with 1082 rows and 2 variables
#' @name SetlowTUV.spct
#'
#' @references
#' \url{http://uv.biospherical.com/Version2/description-Version2-Database3.html} downloaded 2015-02-07
NULL
#' CIE2008 luminous efficiency function (2-deg) (photopic human vision)
#'
#' A dataset containing the wavelengths at a 1 nm interval.
#' Tabulated values for quantum luminous efficiency according to CIE2008 for 2 degrees.
#'
#' The variables are as follows:
#'
#' \itemize{
#'   \item w.length (nm)
#'   \item s.q.response
#' }
#'
#' @docType data
#' @keywords datasets
#' @format A response.spct object with 441 rows and 2 variables
#' @name CIE2008_lef2deg.spct
#' @references
#' \url{http://www.cvrl.org/} downloaded on 2015-01-24
NULL
#' CIE1924 luminous efficiency function (photopic human vision)
#'
#' A dataset containing the wavelengths at a 1 nm interval.
#' Tabulated values for quantum luminous efficiency according to CIE1924.
#'
#' The variables are as follows:
#'
#' \itemize{
#'   \item w.length (nm)
#'   \item s.q.response
#' }
#'
#' @docType data
#' @keywords datasets
#' @format A response.spct object with 471 rows and 2 variables
#' @name CIE1924_lef.spct
#'
#' @note This luminous efficiency function understimates the renponse to short
#'   wavelengths.
#' @references
#'  \url{http://www.cvrl.org/} downloaded on 2015-01-24
NULL
#' Luminous efficiency function (scotopic human vision)
#'
#' A dataset containing the wavelengths at a 1 nm interval. Tabulated values for
#' quantum luminous efficiency at low light levels according to CIE1951.
#'
#' The variables are as follows:
#'
#' \itemize{
#'   \item w.length (nm)
#'   \item s.q.response
#' }
#'
#' @docType data
#' @keywords datasets
#' @format A response.spct object with 401 rows and 2 variables
#' @name CIE1951_scotopic_lef.spct
#' @references
#' \url{http://www.cvrl.org/} downloaded on 2015-01-24
NULL
#' Photopic sensitivity of the human eye
#'
#' Constant value used in the definition of Lumen
#' 1 Lumen is equal to 683 W at 555 nm
#'
#' A single numeric value
#'
#' @format A single numeric value
#' @name photopic_sensitivity
#' @export
#'
photopic_sensitivity <- 683
#' Scotopic sensitivity of the human eye
#'
#' Constant value for human vision under very weak illumination
#' 1 Lumen is equal to 1699 W at 507 nm
#'
#' A single numeric value
#'
#' @format A single numeric value
#' @name scotopic_sensitivity
#' @export
#'
scotopic_sensitivity <- 1699
