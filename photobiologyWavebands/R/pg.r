#' Definition of PG weighted waveband
#'
#' Plant growth BSWF
#'
#' @param norm normalization wavelength (nm)
#' @param w.low short-end boundary wavelength (nm)
#' @param w.high long-end boundary wavelength (nm)
#'
#' @return a waveband object wavelength defining wavelength range, weighting
#'   function and normalization wavelength.
#'
#' @references [1] Flint, S. and Caldwell M. M. (2003) A biological spectral
#' weighting function for ozone depletion research with higher plants
#' Physiologia Plantarum, 2003, 117, 137-144
#'
#' [2] Micheletti, M. I.; Piacentini, R. D. & Madronich, S. (2003) Sensitivity
#' of Biologically Active UV Radiation to Stratospheric Ozone Changes: Effects
#' of Action Spectrum Shape and Wavelength Range Photochemistry and
#' Photobiology, 78, 456-461
#'
#' [3] \url{http://cprm.acd.ucar.edu/Models/TUV/}
#'
#' @note In the original publication [1], the long-end wavelength boundary is
#'   not specified. the longest wavelength at which the plant response was
#'   measured is 366 nm. From the data there is no evidence that action would
#'   immediately drop to zero at longer wavelengths. We have used in earlier
#'   versions the same value as used by the 'NSF Polar Programs UV Monitoring
#'   Network' as described in
#'   \url{http://uv.biospherical.com/Version2/description-Version2-Database3.html}.
#'    Now we keep 390 nm as our default value, but make if possible for the user
#'   to set a different wavelength. To reproduce the output of the TUV
#'   simulation model [3] version 5.0 set \code{w.high = 366}.
#'
#'   In contrast to the NSF Network, for example, the programme TUV uses 366 nm
#'   as the limit, so for comparing results one may need to adjust the value of
#'   this parameter. The effect on the RAF and doses of changing this wavelength
#'   boundary is substantial, as discussed by Micheletti et al. [2].
#'
#' @export
#' @seealso \code{\link{GEN_G}} \code{\link{GEN_T}} \code{\link{GEN_M}} and
#'  \code{\link[photobiology]{waveband}}
#' @examples
#' PG()
#' PG(300)
#'
#' @family BSWF weighted wavebands
#'
PG <- function(norm=300, w.low=275, w.high=390) {
  new_waveband(w.low=w.low, w.high=w.high, weight="SWF",
               SWF.q.fun=PG_q_fun, SWF.norm=299.925658,
               norm=norm,
               wb.name=paste("PG", as.character(norm), sep="."), wb.label="PG")
}
