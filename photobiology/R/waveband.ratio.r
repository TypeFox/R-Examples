#' Photon or energy ratio
#'
#' This function gives the (energy or photon) irradiance ratio between two given
#' wavebands of a radiation spectrum.
#'
#' @param w.length numeric Vector of wavelengths (nm)
#' @param s.irrad numeric Corresponding vector of spectral (energy) irradiances
#'   (W m-2 nm-1)
#' @param w.band.num waveband
#' @param w.band.denom waveband
#' @param unit.out.num character Allowed values "energy", and "photon", or its
#'   alias "quantum"
#' @param unit.out.denom character Allowed values "energy", and "photon", or its
#'   alias "quantum"
#' @param unit.in character Allowed values "energy", and "photon", or its alias
#'   "quantum"
#' @param check.spectrum logical Flag indicating whether to sanity check input
#'   data, default is TRUE
#' @param use.cached.mult logical Flag indicating whether multiplier values
#'   should be cached between calls
#' @param use.hinges logical Flag indicating whether to use hinges to reduce
#'   interpolation errors
#'
#' @return a single numeric value giving the ratio
#'
#' @export
#' @examples
#' # photon:photon ratio
#' with(sun.data,
#'      waveband_ratio(w.length, s.e.irrad,
#'                     new_waveband(400,500),
#'                     new_waveband(400,700), "photon"))
#' # energy:energy ratio
#' with(sun.data,
#'      waveband_ratio(w.length, s.e.irrad,
#'                     new_waveband(400,500),
#'                     new_waveband(400,700), "energy"))
#' # energy:photon ratio
#' with(sun.data,
#'      waveband_ratio(w.length, s.e.irrad,
#'                     new_waveband(400,700),
#'                     new_waveband(400,700),
#'                     "energy", "photon"))
#' # photon:photon ratio waveband : whole spectrum
#' with(sun.data,
#'      waveband_ratio(w.length, s.e.irrad,
#'                     new_waveband(400,500),
#'                     unit.out.num="photon"))
#' # photon:photon ratio of whole spectrum should be equal to 1.0
#' with(sun.data,
#'      waveband_ratio(w.length, s.e.irrad,
#'      unit.out.num="photon"))
#'
waveband_ratio <-
  function(w.length, s.irrad,
           w.band.num = NULL, w.band.denom = NULL,
           unit.out.num = NULL, unit.out.denom = unit.out.num,
           unit.in = "energy",
           check.spectrum = TRUE,
           use.cached.mult = FALSE,
           use.hinges = getOption("photobiology.use.hinges",
                                  default = NULL) ) {
    # We duplicate code from irradiance() here to avoid repeated checks
    # and calculations on the same data
    #
    # what output? seems safer to not have a default here
    if (is.null(unit.out.num) || is.null(unit.out.denom)) {
      warning("'unit.out.num' has no default value")
      return(NA)
    }
    # make code a bit simpler further down
    if (unit.in == "quantum") {unit.in <- "photon"}
    # sanity check for wavelengths
    if (check.spectrum && !check_spectrum(w.length, s.irrad)) {
      return(NA)
    }
    # if the waveband for numerator is undefined then use
    # the whole wavelength range of the spectrum for numerator
    if (is.null(w.band.num)) {
      w.band.num <- new_waveband(min(w.length),max(w.length))
      warning("'w.band.num' not supplied, using whole range of data instead.")
    }
    # if the waveband for denominator is undefined then use
    # the whole wavelength range of the spectrum for denominator
    if (is.null(w.band.denom)) {
      w.band.denom <- new_waveband(min(w.length),max(w.length))
      warning("'w.band.denom' not supplied, using whole range of data instead.")
    }
    # choose whether to use hinges or not
    # if the user has specified its value, we leave it alone
    # but if it was not requested, we decide whether to use
    # it or not based of the wavelength resolution of the
    # spectrum. This will produce small errors for high
    # spectral resulution data, and speed up the calculations
    # a lot in such cases
    if (is.null(use.hinges)) {
      use.hinges <- auto_hinges(w.length)
    }
    # if the w.band.num and/or w.band.denom include 'hinges' we insert them
    # it is o.k. to have hinges unsorted!
    # in new_waveband() NULL hinges are replaced with numeric(0)
    if (use.hinges) {
      merged.hinges <- c(w.band.denom$hinges, w.band.num$hinges)
      if (length(merged.hinges) > 0) {
        new.data <- l_insert_hinges(x = w.length, y = s.irrad, merged.hinges)
        w.length <- new.data$x
        s.irrad <- new.data$y
      }
    }
    # calculate the multipliers
    mult.num <- calc_multipliers(w.length, w.band.num,
                                 unit.out.num, unit.in,
                                 use.cached.mult = use.cached.mult)
    mult.denom <- calc_multipliers(w.length, w.band.denom,
                                   unit.out.denom, unit.in,
                                   use.cached.mult = use.cached.mult)

    # calculate weighted spectral irradiance
    irrad.num <- integrate_xy(w.length, s.irrad * mult.num)
    irrad.denom <- integrate_xy(w.length, s.irrad * mult.denom)

    return(irrad.num / irrad.denom)
  }
