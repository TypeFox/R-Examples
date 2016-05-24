#' Convert spectral photon irradiance into spectral energy irradiance
#'
#' For example an spectrum [mol s-1 m-2 nm-1] is converted into a spectrum [W
#' m-2 nm-1]
#'
#' @param w.length numeric Vector of wavelengths (nm)
#' @param s.qmol.irrad numeric Corresponding vector of spectral photon
#'   irradiances
#'
#' @return A numeric vector of spectral (energy) irradiances
#' @export
#'
#' @examples
#' with(sun.spct, as_energy(w.length, s.q.irrad))
#'
#' @family quantity conversion functions
#'
as_energy <- function(w.length, s.qmol.irrad){
  return(s.qmol.irrad / e2qmol_multipliers(w.length))
}
