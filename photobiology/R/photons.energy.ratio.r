#' Photon:energy ratio
#'
#' This function gives the photons:energy ratio between for one given waveband
#' of a radiation spectrum.
#'
#' @param w.length numeric array of wavelength (nm)
#' @param s.irrad numeric array of spectral (energy) irradiances (W m-2 nm-1)
#' @param w.band waveband
#' @param unit.in character Allowed values "energy", and "photon", or its alias
#'   "quantum"
#' @param check.spectrum logical Flag telling whether to sanity check input
#'   data, default is TRUE
#' @param use.cached.mult logical Flag telling whether multiplier values should
#'   be cached between calls
#' @param use.hinges logical Flag telling whether to use hinges to reduce
#'   interpolation errors
#'
#' @return A single numeric value giving the ratio moles-photons per Joule.
#'
#' @export
#' @examples
#' # photons:energy ratio
#' with(sun.data, photons_energy_ratio(w.length, s.e.irrad, new_waveband(400,500)))
#' # photons:energy ratio for whole spectrum
#' with(sun.data, photons_energy_ratio(w.length, s.e.irrad))
#'
#' @family photon and energy ratio functions
#'
photons_energy_ratio <- function(w.length, s.irrad,
                           w.band=NULL,
                           unit.in="energy",
                           check.spectrum=TRUE,
                           use.cached.mult = FALSE,
                           use.hinges = getOption("photobiology.use.hinges", default=NULL) ){
  return(waveband_ratio(w.length, s.irrad, w.band, w.band,
                        unit.out.num="photon", unit.out.denom="energy",
                        unit.in=unit.in,
                        check.spectrum=check.spectrum,
                        use.cached.mult=use.cached.mult,
                        use.hinges = use.hinges))
 }
