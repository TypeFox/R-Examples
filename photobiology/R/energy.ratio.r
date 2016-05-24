#' Energy:energy ratio
#'
#' This function gives the energy ratio between two given wavebands of a
#' radiation spectrum.
#'
#' @param w.length numeric Vector of wavelengths (nm)
#' @param s.irrad numeric Corresponding of spectral (energy) irradiances (W m-2
#'   nm-1)
#' @param w.band.num waveband
#' @param w.band.denom waveband
#' @param unit.in character Allowed values "energy", and "photon", or its alias
#'   "quantum"
#' @param check.spectrum logical Flag indicating whether to sanity check input
#'   data, default is TRUE
#' @param use.cached.mult logical Flag indicating whether multiplier values
#'   should be cached between calls
#' @param use.hinges logical Flag indicating whether to use hinges to reduce
#'   interpolation errors
#'
#' @return a single numeric value giving the unitless ratio
#'
#' @export
#' @examples
#' # energy:energy ratio
#' with(sun.data,
#'      energy_ratio(w.length, s.e.irrad, new_waveband(400,500), new_waveband(400,700)))
#' # energy:energy ratio waveband : whole spectrum
#' with(sun.data, energy_ratio(w.length, s.e.irrad, new_waveband(400,500)))
#' # energy:energy ratio of whole spectrum should be equal to 1.0
#' with(sun.data, energy_ratio(w.length, s.e.irrad))
#'
#' @family photon and energy ratio functions
#'
energy_ratio <- function(w.length, s.irrad,
                         w.band.num = NULL, w.band.denom = NULL,
                         unit.in = "energy",
                         check.spectrum = TRUE,
                         use.cached.mult = FALSE,
                         use.hinges = NULL) {
  return(waveband_ratio(w.length, s.irrad, w.band.num, w.band.denom,
                        unit.out.num = "energy", unit.out.denom = "energy",
                        unit.in = unit.in,
                        check.spectrum = check.spectrum,
                        use.cached.mult = use.cached.mult,
                        use.hinges = use.hinges))
}
