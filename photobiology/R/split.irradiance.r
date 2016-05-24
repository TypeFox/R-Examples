#' Energy or photon irradiance for split spectrum regions
#'
#' This function returns the energy or photon irradiance for a series of
#' contiguous wavebands from a radiation spectrum. The returned values can be
#' either absolute or relative to their sum.
#'
#' @param w.length numeric Vector of wavelengths (nm)
#' @param s.irrad numeric Corresponding vector of spectral (energy) irradiances (W m-2 nm-1)
#' @param cut.w.length numeric Vector of wavelengths (nm)
#' @param unit.out character Allowed values "energy", and "photon",
#'   or its alias "quantum"
#' @param unit.in character Allowed values "energy", and "photon",
#'   or its alias "quantum"
#' @param scale a character A string indicating the scale used for the returned
#'   values ("absolute", "relative", "percent")
#' @param check.spectrum logical Flag indicating whether to sanity check input data,
#'   default is TRUE
#' @param use.cached.mult logical Flag indicating whether multiplier values should be
#'   cached between calls
#' @param use.hinges logical Flag indicating whether to use hinges to reduce
#'   interpolation errors
#'
#' @return a numeric array of irradiances with no change in scale factor: [W m-2
#'   nm-1] -> [mol s-1 m-2] or relative values (fraction of one) if scale =
#'   "relative" or scale = "percent"
#'
#' @export
#' @examples
#' with(sun.data,
#'       split_irradiance(w.length, s.e.irrad,
#'                        cut.w.length = c(300, 400, 500, 600, 700),
#'                        unit.out = "photon"))
#' @note The last three parameters control speed optimizations. The defaults
#'   should be suitable in mosts cases. If you set \code{check.spectrum=FALSE}
#'   then you should call \code{\link{check_spectrum}} at least once for your
#'   spectrum before using any of the other functions. If you will use
#'   repeatedly the same SWFs on many spectra measured at exactly the same
#'   wavelengths you may obtain some speed up by setting
#'   \code{use.cached.mult=TRUE}. However, be aware that you are responsible for
#'   ensuring that the wavelengths are the same in each call, as the only test
#'   done is for the length of the \code{w.length} vector.
#'
#' @family split a spectrum into regions functions
#'
split_irradiance <- function(w.length, s.irrad,
                             cut.w.length=range(w.length),
                             unit.out=getOption("photobiology.base.unit", default="energy"),
                             unit.in="energy",
                             scale="absolute",
                             check.spectrum=TRUE,
                             use.cached.mult = FALSE,
                             use.hinges=getOption("photobiology.use.hinges", default=NULL) )
{
  # what output? seems safer to not have a default here
  if (is.null(unit.out)){
    warning("'unit.out' has no default value")
    return(NA)
  }
  # make code a bit simpler further down
  if (unit.in=="quantum") {unit.in <- "photon"}
  # sanity check for spectral data and abort if check fails
  if (check.spectrum && !check_spectrum(w.length, s.irrad)) {
    return(NA)
  }
  # clean the cut point wavelengths
  if (is.null(cut.w.length)){
    cut.w.length=range(w.length)
  } else {
    cut.w.length <- unique(sort(cut.w.length))
    len.cut.w.length <- length(cut.w.length)
    if (len.cut.w.length < 2) {
      warning("End cut points have no default, but only one cut point supplied")
      return(NA)
    } else if (cut.w.length[1] > max(w.length) || cut.w.length[len.cut.w.length] < min(w.length)) {
      warning("All cut points are outside the range of the spectral data")
      return(NA)
    }
    if (cut.w.length[1] < min(w.length)) {
      warning("Shortest cut point(s) outside spectral data range")
      i <- 2
      while(i < len.cut.w.length && cut.w.length[i] < min(w.length)) i <- i + 1
      cut.w.length <- c(min(w.length), cut.w.length[i:len.cut.w.length])
      len.cut.w.length <- length(cut.w.length)
    }
    if (cut.w.length[len.cut.w.length] > max(w.length)) {
      warning("Longest cut point(s) outside spectral data range")
      j <- len.cut.w.length - 1
      while(j > 2 && cut.w.length[j] > max(w.length)) j <- j - 1
      cut.w.length <- c(cut.w.length[1:j], max(w.length))
      len.cut.w.length <- length(cut.w.length)
    }
    if (len.cut.w.length < 2) {
      stop("Failed assertion in split.irradiance. You've found a bug!")
    }
  }
  w.band.num <- length(cut.w.length) - 1
  w.bands <- list(new_waveband(cut.w.length[1], cut.w.length[2]))
  if (w.band.num > 1) for (i in 2:w.band.num) {
    w.bands <- c(w.bands, list(new_waveband(cut.w.length[i], cut.w.length[i+1])))
  }
  irrads <- irradiance(w.length, s.irrad, w.bands, unit.out=unit.out, unit.in=unit.in,
                       check.spectrum=FALSE, use.cached.mult=use.cached.mult, use.hinges=use.hinges)
  if (scale == "relative") {
    irrads <- irrads / sum(irrads)
  } else if (scale == "percent") {
    irrads <- irrads / sum(irrads) * 100
  }
  return(irrads)
}
