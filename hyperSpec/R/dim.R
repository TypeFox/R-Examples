##' The Number of Rows (Spectra), Columns, and Data Points per Spectrum of an
##' hyperSpec Object)
##'
##' \code{ncol} returns the number of columns in \code{x@@data}. I.e. the number
##' of columns with additional information to each spectrum (e.g. "x", "y",
##' \dots{}) + 1 (for column \code{spc} containing the spectra).
##' @rdname dim
##' @docType methods
##' @param x a \code{hyperSpec} object
##' @author C. Beleites
##' @seealso \code{\link[base]{ncol}}
##' 
##' @return \code{nrow}, \code{ncol}, \code{nwl}, and \code{length}, return an
##'   \code{integer}.
##' @export 
##' @examples
##' 
##' ncol (chondro)
setMethod ("ncol", signature = signature ("hyperSpec"), function (x){
  validObject (x)

  ncol (x@data)
})

##'
##' \code{nrow} yields the number of rows in \code{x@@data}, i.e. the number of
##' spectra in the \code{hyperSpec} object.
##' 
##' @rdname dim
##' @seealso \code{\link[base]{nrow}}
##' @export
##' @examples
##' nrow (chondro)
setMethod ("nrow", signature = signature ("hyperSpec"), function (x){
  validObject (x)

  nrow (x@data)
})

##'
##' \code{nwl} returns the number of columns in \code{x@@data$spc}, i.e. the
##' length of each spectrum.
##' 
##' @rdname dim
##' @aliases nwl
##' @export
##' @examples
##' 
##' nwl  (chondro)
nwl <- function (x){
  chk.hy (x)
  validObject (x)
  
  ncol (x@data$spc)
}



##'
##' \code{dim} returns all three values in a vector.
##' 
##' 
##' @rdname dim
##' @return
##' 
##' \code{dim} returns a vector of length 3.
##' @seealso \code{\link[base]{dim}}
##' @keywords methods
##' @export 
##' @examples
##' dim (chondro)
setMethod ("dim", signature = signature ("hyperSpec"), function (x){
  validObject (x)
  c (nrow = nrow (x@data), ncol = ncol (x@data), nwl = ncol (x@data$spc))
})

##' 
##' \code{length} is a synonym for \code{nrow}. It is supplied so that
##' \code{seq_along (x)} returns a sequence to index each spectrum.
##' @rdname dim
##' @seealso \code{\link[base]{length}}
##' @export 
##' @examples
##' length (chondro)
setMethod ("length", signature = signature ("hyperSpec"), function (x) {
  validObject (x)
  nrow (x@data)
})


