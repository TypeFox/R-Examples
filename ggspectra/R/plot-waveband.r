#' Plot a waveband as a spectrum.
#'
#' This function returns a ggplot object with an annotated plot of a
#' waveband object.
#'
#' @note Note that scales are expanded so as to make space for the annotations.
#'   The object returned is a ggplot object, and can be further manipulated.
#'
#' @param x a waveband object
#' @param ... other arguments passed to plot.response_spct()
#' @param w.length numeric vector of wavelengths (nm)
#' @param range an R object on which range() returns a vector of length 2, with
#'   min annd max wavelengths (nm)
#' @param fill value to use as response for wavelngths outside the waveband
#'   range
#' @param unit.in the type of unit we assume as reference "energy" or "photon"
#'   based
#' @param annotations a character vector
#' @param wb.trim logical
#' @param norm numeric normalization wavelength (nm) or character string "max"
#'   for normalization at the wavelength of highest peak.
#'
#' @return a \code{ggplot} object.
#'
#' @note Effectiveness spectra are plotted expressing the spectral effectiveness
#' either as $1 mol^{-1} nm$ photons of $1 J^{-1} nm$ which can selected through
#' formal argument \code{unit.out}. The value of \code{unit.in} has no effect on
#' the result when uisng BSWFs, as BSWFs are defined based on a certain base of
#' expression, which is enforced. In contrast, for wavebands which only define a
#' wavelength range, changing the assumed reference irradiance, changes the
#' responsivity according to Plank's law.
#'
#' This function creates a response_spct object from the waveband object and
#' plots it. Unused arguments are passed along, which means that other plot
#' aspects can be controlled by providing arguments for the plot method of the
#' response_spct class.
#'
#' @keywords hplot
#'
#' @export
#'
#' @family plot functions
#'
#' @examples
#' library(photobiology)
#' plot(waveband(c(400, 500)))
#'
plot.waveband <-
  function(x,
           ...,
           w.length = NULL,
           range = c(280, 800),
           fill = 0,
           unit.in = getOption("photobiology.radiation.unit", default = "energy"),
           annotations = getOption("photobiology.plot.annotations",
                                   default = c("colour.guide", "boxes", "labels")),
           wb.trim = TRUE,
           norm = NULL) {
  w.band <- x
  if (!is.waveband(w.band)) {
    return(ggplot())
  }
  if (!is.null(w.length)) {
    w.length <- unique(sort(w.length, na.last = NA))
  }
  if (is.null(range)) {
    if (is.null(w.length) || length(w.length) < 2) {
      range <- range(w.band)
    }
  } else {
    range <- range(range)
  }
  w.length <- w.length[w.length > range[1] & w.length < range[2]]
  if (is.null(w.length)) {
    w.length <- seq(range[1], range[2], length.out = 200)
  } else if (length(w.length) < 200) {
    range <- range(w.length)
    w.length <- seq(range[1], range[2], length.out = 200)
  }
  if (!is.null(w.band$hinges) & length(w.band$hinges)>0) {
    hinges <- with(w.band, hinges[hinges > range[1] & hinges < range[2]])
    w.length <- c(w.length, hinges)
  }
  w.length <- unique(sort(w.length))
  s.response <-
    calc_multipliers(w.length, w.band,
                     unit.out=unit.in, unit.in=unit.in,
                     use.cached.mult=getOption("photobiology.use.cached.mult",
                                               default = FALSE), fill = fill)
  if (is.null(norm)) {
    if (!is.null(w.band$norm)) {
      norm <- w.band$norm
    } else {
      norm <- "max"
    }
  }
  if (unit.in == "energy") {
    spct <- response_spct(w.length = w.length, s.e.response = s.response)
  } else if (unit.in %in% c("photon", "quantum")) {
    spct <- response_spct(w.length = w.length, s.q.response = s.response)
  }
  if (is_effective(w.band)) {
    w.band.range <-
      waveband(w.band,
               wb.name = paste("Range of", labels(w.band)[["label"]]))
  } else {
    w.band.range <- w.band
  }
  out.ggplot <- plot(spct, w.band=w.band.range, annotations = annotations,
                     wb.trim = wb.trim, norm = norm, ...)
  if ("title" %in% annotations) {
    out.ggplot <- out.ggplot + labs(title = deparse(substitute(x)))
  }
  return(out.ggplot)
}