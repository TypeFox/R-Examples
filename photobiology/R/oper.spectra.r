#' Binary operation on two spectra, even if the wavelengths values differ
#'
#' The wavelength vectors of the two spectra are merged, and the missing
#' spectral values are calculated by interpolation. After this, the two spectral
#' values at each wavelength are added.
#'
#' @param w.length1 numeric array of wavelength (nm)
#' @param w.length2 numeric array of wavelength (nm)
#' @param s.irrad1 a numeric array of spectral values
#' @param s.irrad2 a numeric array of spectral values
#' @param trim a character string with value "union" or "intersection"
#' @param na.rm a logical value, if TRUE, not the default, NAs in the input are
#'   replaced with zeros
#' @param bin.oper a function defining a binary operator (for the usual math
#'   operators enclose argument in backticks)
#' @param ... additional arguments (by name) passed to bin.oper
#' @return a dataframe with two numeric variables \item{w.length}{A numeric
#'   vector with the wavelengths (nm) obtained by "fusing" w.length1 and
#'   w.length2. w.length contains all the unique vales, sorted in ascending
#'   order.} \item{s.irrad}{A numeric vector with the sum of the two spectral
#'   values at each wavelength.}
#' @details If trim=="union" spectral values are calculated for the whole range
#'   of wavelengths covered by at least one of the input spectra, and missing
#'   values are set in each input spectrum to zero before addition. If
#'   trim=="intersection" then the range of wavelengths covered by both input
#'   spectra is returned, and the non-overlaping regions discarded. If
#'   w.length2==NULL, it is assumed that both spectra are measured at the same
#'   wavelengths, and a simple addition is used, ensuring fast calculation.
#' @export
#'
#' @examples
#' 
#' head(sun.data)
#' result.data <-
#'     with(sun.data,
#'          oper_spectra(w.length, w.length, s.e.irrad, s.e.irrad, bin.oper=`+`))
#' head(result.data)
#' tail(result.data)
#' my_fun <- function(e1, e2, k) {return((e1 + e2) / k)}
#' result.data <-
#'     with(sun.data,
#'         oper_spectra(w.length, w.length, s.e.irrad, s.e.irrad, bin.oper=my_fun, k=2))
#' head(result.data)
#' tail(result.data)
#'
oper_spectra <- function(w.length1, w.length2=NULL, s.irrad1, s.irrad2, trim="union", na.rm=FALSE, bin.oper=NULL, ...) {
  if (na.rm) {
    ifelse(!is.na(s.irrad1), s.irrad1, 0.0)
    ifelse(!is.na(s.irrad2), s.irrad2, 0.0)
  }
  if (is.null(w.length2)) {
    if (length(s.irrad1) == length(s.irrad2) & length(w.length1) == length(s.irrad1)){
      s.irrad.result <- bin.oper(s.irrad1, s.irrad2, ...)
      w.length <- w.length1
    } else {
      stop("Mismatch in the length of input vectors")
    }
    invisible(dplyr::data_frame(w.length, s.irrad=s.irrad.result))
  }
  if (length(w.length2) != length(s.irrad2) | length(w.length1) != length(s.irrad1)){
    stop("Mismatch in the length of input vectors")
  }
  else if (trim=="union") {
    wl.low <- min(w.length1[1],w.length2[1])
    wl.hi <- max(w.length1[length(w.length1)],w.length2[length(w.length2)])
  }
  else if (trim=="intersection") {
    wl.low <- max(w.length1[1],w.length2[1])
    wl.hi <- min(w.length1[length(w.length1)],w.length2[length(w.length2)])
  }
  else {
    warning("illegal value for 'trim' argument")
    invisible(NA)
  }
  w.length <- c(w.length1[w.length1 >= wl.low & w.length1 <= wl.hi], w.length2[w.length2 >= wl.low & w.length2 <= wl.hi])
  w.length <- sort(w.length)
  w.length <- unique(w.length)
  s.irrad1.int <- rep(NA, length(w.length))
  s.irrad2.int <- rep(NA, length(w.length))
  s.irrad1.int <- interpolate_spectrum(w.length1, s.irrad1, w.length, fill=0.0)
  s.irrad2.int <- interpolate_spectrum(w.length2, s.irrad2, w.length, fill=0.0)
  s.irrad.result <- bin.oper(s.irrad1.int, s.irrad2.int, ...)
  invisible(dplyr::data_frame(w.length, s.irrad=s.irrad.result))
}
