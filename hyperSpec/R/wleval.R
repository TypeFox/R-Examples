##' Evaluate function on wavelengths of hyperSpec object
##'
##' This is useful for generating certain types of baseline "reference spectra".
##'
##' @param x hyperSpec object
##' @param ... hyperSpec method: expressions to be evaluated
##' @param normalize.wl function to transorm the wavelengths before evaluating the polynomial (or
##' other function). Use \code{\link[hyperSpec]{normalize01}} to map the wavelength range to the interval [0, 1].
##' @return hyperSpec object containing one spectrum for each expression 
##' @export
##' @seealso \code{\link[hyperSpec]{vanderMonde}} for  polynomials,
##'
##' \code{\link[hyperSpec]{normalize01}} to normalize the wavnumbers before evaluating the function
##' @author C. Beleites
##' @examples
##' plot (wl.eval (laser, exp = function (x) exp (-x)))
wl.eval <- function (x, ..., normalize.wl = I){
  chk.hy (x)
  validObject (x)

  fun <- list (...)

  wl <- normalize.wl (x@wavelength)
  
  x <- decomposition (x, t (sapply (fun, function (f) f (wl))), scores = FALSE)
  x$.f <- if (is.null (names (fun)))
              rep (NA, length (fun))
          else
              names (fun)
  x
}

.test (wl.eval) <- function (){
  x <- runif (10, min = -1e3, max = 1e3)
  
  checkEqualsNumeric (min (normalize01 (x)), 0)
  checkEqualsNumeric (max (normalize01 (x)), 1)

  checkEqualsNumeric (normalize01 (x), (x - min (x)) / diff (range (x)))
}

