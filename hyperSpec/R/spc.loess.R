##' loess smoothing interpolation for spectra
##' Spectra can be smoothed and interpolated on a new wavelength axis using
##' \code{\link[stats]{loess}}.
##' 
##' Applying \code{\link[stats]{loess}} to each of the spectra, an interpolation onto a new
##' wavelength axis is performed.  At the same time, the specta are smoothed in order to increase the
##' signal : noise ratio. See \code{\link[stats]{loess}} and \code{\link[stats]{loess.control}} on
##' the parameters that control the amount of smoothing.
##' 
##' @param spc the \code{hyperSpec} object
##' @param newx wavelengh axis to interpolate on
##' @param enp.target,surface,... parameters for \code{\link[stats]{loess}} and
##' \code{\link[stats]{loess.control}}. 
##' @return a new \code{hyperspec} object.
##' @rdname spc-loess
##' @export
##' @author C. Beleites
##' @seealso \code{\link[stats]{loess}}, \code{\link[stats]{loess.control}}
##' @keywords manip datagen
##' @examples
##' 
##' plot (flu, col = "darkgray")
##' plot (spc.loess(flu, seq (420, 470, 5)), add = TRUE, col = "red")
##' 
##' flu [[3, ]] <- NA_real_
##' smooth <- spc.loess(flu, seq (420, 470, 5))
##' smooth [[, ]]
##' plot (smooth, add = TRUE, col = "blue")
##' 
spc.loess <- function (spc, newx, enp.target = nwl (spc) / 4,
                       surface = "direct", ...){

  .loess <- function (y, x){
    if (all (is.na (y)))
      NA
    else
      loess (y ~ x, enp.target = enp.target, surface = surface, ...)
  }

  .predict <-  function (loess, x){
    if (!is (loess, "loess") && is.na (loess))
      rep (NA_real_, length (x))
    else
      predict (loess, x)
  }
  
  chk.hy (spc)
  validObject (spc)

  loess <- apply (t (spc[[]]), 2, .loess, spc@wavelength)

  spc@data$spc <- t (sapply (loess, .predict, newx))
  .wl(spc) <- newx

  if (any (is.na (spc@data$spc)))
    warning ("NAs were generated. Probably newx was outside the spectral range covered by spc.")

  spc
}
