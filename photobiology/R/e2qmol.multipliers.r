#' Calculate energy to quantum (mol) multipliers
#'
#' Gives multipliers as a function of wavelength, for converting from energy to
#' photon (quantum) molar units.
#'
#' @param w.length numeric Vector of wavelengths (nm)
#'
#' @return A numeric array of multipliers
#'
#' @export
#' @examples
#' with(sun.data, e2qmol_multipliers(w.length))
#'
#' @family quantity conversion functions
#'
e2qmol_multipliers <- function(w.length){
  return(e2quantum_multipliers(w.length, molar=TRUE))
}
