##' Smoothing splines
##'
##' Spectral smoothing by splines
##' @title Spectral smoothing by splines
##' @param spc hyperSpec object
##' @param newx  wavelengh axis to interpolate on
##' @param ... further parameters handed to \code{\link[stats]{smooth.spline}}
##' @return hyperSpec object containing smoothed spectra
##' @rdname spc-spline
##' @author Claudia Beleites
##' @seealso \code{\link[hyperSpec]{spc.loess}}
##' 
##' \code{\link[stats]{smooth.spline}}
##' @note This function is still experimental
##' @export
##' @examples
##' p <- paracetamol [,,2200 ~ max]
##' plot (p, col = "gray")
##' smooth <- spc.smooth.spline (p [,, c (2200 ~ 2400, 2500 ~ 2825, 3150 ~ max)],
##'                              wl (paracetamol [,, 2200 ~ max]),
##'                              df = 4, spar = 1)
##' plot (smooth, col = "red", add = TRUE)
##' 
##' plot (p - smooth)
##' 
spc.smooth.spline <- function (spc, newx = wl (spc), ...){

  .spline <- function (x, y, newx){
    pts <- ! is.na (y)
    fit <- smooth.spline (x [pts], y [pts], ...)$fit
    predict (fit, newx, deriv = 0)$y
  }

  spc <- orderwl (spc) #includes chk.hy and validObject

  newspc <- matrix (NA_real_, ncol = length (newx), nrow = nrow (spc))
  i <- rowSums (is.na (spc@data$spc)) < nwl (spc)

  newspc [i, ] <- t (apply (spc@data$spc [i,, drop = FALSE], 1,
                            .spline, x = spc@wavelength, newx = newx))

  if (any (is.na (newspc [i, ])))
      warning ("NAs generated. Probably newx was outside the spectral range covered by spc.")

  spc@data$spc <- newspc
  .wl(spc) <- newx
  
  validObject (spc)
  

  spc
}

