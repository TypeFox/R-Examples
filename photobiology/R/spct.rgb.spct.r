#' RGB color values
#'
#' This function returns the RGB values for a source spectrum.
#'
#' @param spct an object of class "source_spct"
#' @param sens a chroma_spct object with variables w.length, x, y, and z, giving
#'   the CC or CMF definition (default is the proposed human CMF according to
#'   CIE 2006.)
#' @param color.name character string for naming the rgb color definition
#'
#' @return A color defined using \code{rgb()}. The numeric values of the RGB
#'   components can be obtained
#'
#' @export
#' @examples
#' rgb_spct(sun.spct)
#'
#' @family color functions
#'
rgb_spct <-
  function(spct, sens=photobiology::ciexyzCMF2.spct, color.name=NULL){
    if (is.source_spct(spct)) {
      return(s_e_irrad2rgb(spct$w.length, spct$s.e.irrad, sens=sens, color.name=color.name))
    } else {
      return(NA)
    }
  }
