##' Function evaluation on hyperSpec objects
##'
##' vandermonde generates van der Monde matrices, the hyperSpec method generates a hyperSpec object
##' containing the van der Monde matrix of the wavelengths of a hyperSpec object.
##'
##' It is often numerically preferrable to map \code{wl (x)} to [0, 1], see the example.
##' 
##' @param x object to evaluate the polynomial on 
##' @param order of the polynomial
##' @rdname vanderMonde
##' @return van der Monde matrix
##' @author C. Beleites
##' @export
##' @include hyperspec-package.R
vanderMonde <- function (x, order, ...){
  if (nargs () > 2)
    stop ('Unknown arguments: ', names (c (...)))
    
  outer (x, 0 : order, `^`)
}

##' @noRd
setGeneric ("vanderMonde")

##' @param normalize.wl function to transorm the wavelengths before evaluating the polynomial (or
##' other function). \code{\link[hyperSpec]{normalize01}} maps the wavelength range to the interval
##' [0, 1]. Use \code{\link[base]{I}} to turn off. 
##' @param ... hyperSpec method: further arguments to \code{\link{decomposition}}
##' @return hyperSpec method: hyperSpec object containing van der Monde matrix as spectra and an additional column ".vdm.order" giving the order of each spectrum (term).
##' @rdname vanderMonde
##' @seealso \code{\link[hyperSpec]{wl.eval}} for calculating arbitrary functions of the wavelength,
##'
##' \code{\link[hyperSpec]{normalize01}} 
##' @export
##' @examples
##' plot (vanderMonde (flu, 2))
##' plot (vanderMonde (flu, 2, normalize.wl = I))
##'
##' 
setMethod ("vanderMonde", signature = signature (x = "hyperSpec"),
           function (x, order, ..., normalize.wl = normalize01){
  validObject (x)

  wl <- normalize.wl (x@wavelength)
  
  x <- decomposition (x, t (vanderMonde (wl, order)), scores = FALSE, ...)
  x$.vdm.order <- 0 : order
  x
})


.test (vanderMonde) <- function (){
  checkEqualsNumeric (vanderMonde (c (1 : 3, 5), 2),
                      matrix (c (1, 1, 1, 1, 1, 2, 3, 5, 1, 4, 9, 25), nrow = 4)
                      )
  checkException (vanderMonde (1, 0, normalize.wl = normalize01))

  checkTrue (chk.hy (vanderMonde (flu, 0)))
  checkTrue (validObject (vanderMonde (flu, 0)))
  
  checkEqualsNumeric (vanderMonde (paracetamol, 3, normalize.wl = I)[[]],
                      t (vanderMonde (wl (paracetamol), 3)))

  checkEqualsNumeric (vanderMonde (paracetamol, 3, normalize.wl = normalize01)[[]],
                      t (vanderMonde (normalize01 (wl (paracetamol)), 3)))
}



