#' Photo:photon ratio
#'
#' This function gives the photon ratio between two given wavebands of a
#' radiation spectrum.
#'
#' @param w.length numeric array of wavelength (nm)
#' @param s.irrad numeric array of spectral (energy) irradiances (W m-2 nm-1)
#' @param w.band.num waveband
#' @param w.band.denom waveband
#' @param unit.in character Allowed values "energy", and "photon",
#'   or its alias "quantum"
#' @param check.spectrum logical Flag telling whether to sanity check input data,
#'   default is TRUE
#' @param use.cached.mult logical Flag telling whether multiplier values should be
#'   cached between calls
#' @param use.hinges logical Flag telling whether to use hinges to reduce
#'   interpolation errors
#'
#' @return a single numeric value giving the unitless ratio
#'
#' @export
#' @examples
#' # photon:photon ratio
#' with(sun.data,
#'      photon_ratio(w.length, s.e.irrad, new_waveband(400,500), new_waveband(400,700)))
#' # photon:photon ratio waveband : whole spectrum
#' with(sun.data, photon_ratio(w.length, s.e.irrad, new_waveband(400,500)))
#' # photon:photon ratio of whole spectrum should be equal to 1.0
#' with(sun.data, photon_ratio(w.length, s.e.irrad))
#'
#' @family photon and energy ratio functions
#'
photon_ratio <- function(w.length, s.irrad,
                         w.band.num=NULL, w.band.denom=NULL,
                         unit.in="energy",
                         check.spectrum=TRUE,
                         use.cached.mult = FALSE,
                         use.hinges=getOption("photobiology.use.hinges", default=NULL) ) {
  return(waveband_ratio(w.length, s.irrad, w.band.num, w.band.denom,
                        unit.out.num="photon", unit.out.denom="photon",
                        unit.in=unit.in,
                        check.spectrum=check.spectrum,
                        use.cached.mult=use.cached.mult,
                        use.hinges=use.hinges))
}
