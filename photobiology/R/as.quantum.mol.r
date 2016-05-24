#' Convert spectral energy irradiance into spectral photon irradiance
#'
#' For example an spectrum [W m-2 nm-1] is converted into a spectrum [mol s-1
#' m-2 nm-1]
#'
#' @param w.length numeric Vector of wavelengths (nm)
#' @param s.e.irrad numeric Corresponding vector of spectral (energy)
#'   irradiances
#'
#' @return a numeric array of spectral photon irradiances
#' @export
#'
#' @examples
#' with(sun.data, as_quantum_mol(w.length, s.e.irrad))
#'
#' @family quantity conversion functions
#'
as_quantum_mol <- function(w.length, s.e.irrad){
  return(s.e.irrad * e2qmol_multipliers(w.length))
}
