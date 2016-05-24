#' Photon irradiance
#'
#' This function returns the photon irradiance for a given waveband of a
#' radiation spectrum, optionally applies a BSWF.
#'
#' @param w.length numeric array of wavelength (nm)
#' @param s.irrad numeric array of spectral irradiances, by default as energy (W
#'   m-2 nm-1)
#' @param w.band waveband
#' @param unit.in character Values recognized "photon" or "energy"
#' @param check.spectrum logical Flag telling whether to sanity check input data,
#'   default is TRUE
#' @param use.cached.mult logical Flag telling whether multiplier values should be
#'   cached between calls
#' @param use.hinges logical Flag telling whether to use hinges to reduce
#'   interpolation errors
#'
#' @return A single numeric value with no change in scale factor: [W m-2 nm-1]
#'   -> [W m-2]
#'
#' @export
#' @examples
#' with(sun.data, photon_irradiance(w.length, s.e.irrad))
#' with(sun.data, photon_irradiance(w.length, s.e.irrad, new_waveband(400,700)))
#'
#' @family irradiance functions
#'
photon_irradiance <-
  function(w.length, s.irrad, w.band=NULL, unit.in="energy",
           check.spectrum=TRUE,
           use.cached.mult = FALSE,
           use.hinges=getOption("photobiology.use.hinges", default=NULL) ) {
    return(irradiance(w.length=w.length, s.irrad=s.irrad, w.band=w.band,
                      unit.out="photon", unit.in=unit.in,
                      check.spectrum=check.spectrum, use.cached.mult=use.cached.mult,
                      use.hinges=use.hinges))
  }
