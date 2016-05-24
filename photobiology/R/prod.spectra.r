#' Multiply two spectra, even if the wavelengths values differ
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
#'
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
#' square.sun.data <-
#'   with(sun.data, prod_spectra(w.length, w.length, s.e.irrad, s.e.irrad))
#' head(square.sun.data)
#' tail(square.sun.data)
#'
prod_spectra <- function(w.length1, w.length2=NULL, s.irrad1, s.irrad2, trim="union", na.rm=FALSE) {
  return(oper_spectra(w.length1=w.length1, w.length2=w.length2,
                      s.irrad1=s.irrad1, s.irrad2=s.irrad2,
                      trim="union", na.rm=FALSE,
                      bin.oper=`*`))
}
