##' Rubberband baseline
##'
##' Baseline with support points determined from a convex hull of the spectrum.
##' @title Rubberband baseline correction
##' @param spc hyperSpec object
##' @param ... further parameters handed to \code{\link[stats]{smooth.spline}}
##' @param upper logical indicating whether the lower or upper part of the hull should be used
##' @param noise noise level to be taken into account
##' @param spline logical indicating whether the baseline should be an interpolating spline through
##' the support points or piecewise linear.
##' @return hyperSpec object containing the baselines
##' @rdname spc-rubberband
##' @author Claudia Beleites
##' @seealso \code{\link[hyperSpec]{spc.fit.poly}}, \code{\link[hyperSpec]{spc.fit.poly.below}}
##'
##' \code{vignette ("baseline")}
##' @note This function is still experimental
##' @export
##' @examples
##' plot (paracetamol [,, 175 ~ 1800])
##' bl <- spc.rubberband (paracetamol [,, 175 ~ 1800], noise = 300, df = 20)
##' plot (bl, add = TRUE, col = 2)
##' plot (paracetamol [,, 175 ~ 1800] - bl)

spc.rubberband <- function (spc, ..., upper = FALSE, noise = 0, spline = TRUE){
  spc <- orderwl (spc)

  if (upper) spc@data$spc <- -spc@data$spc
  
  spc@data$spc <- .rubberband (spc@wavelength, spc@data$spc, 
                               noise = noise, spline = spline, ...)

  if (upper) spc@data$spc <- -spc@data$spc

  spc
}

.rubberband <- function (x, y, noise, spline, ...){
  for (s in seq_len (nrow (y))){
    pts <- chull (x, y [s,])
    imax <- which.max (pts) - 1
    pts <- c (pts [- seq_len (imax)], pts [seq_len (imax)])

    neg <- which (diff (pts) < 0)
    pts <- c (ncol (y), pts [c (1, neg + 1)])
    pts <- sort (unique (pts))

    tmp <- approx (x = x [pts], y = y [s, pts], xout= x, method="linear")$y
    
    if (spline){
      pts <- which (y [s,] <= tmp + noise)

      if (length (pts) > 3)
        tmp <- predict (smooth.spline (x[pts], y[s, pts], ...)$fit, x, 0)$y 
      else 
        tmp <- spline (x [pts], y [s, pts], xout = x)$y
        
    }
    
    y [s, ] <- tmp
    
  }
  
  y  
}
