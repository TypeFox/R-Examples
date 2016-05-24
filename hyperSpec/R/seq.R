
##' Sequence generation along spectra or wavelengths
##' This function generates sequences along the spectra (rows) or wavelengths of hyperSpec objects.
##' 
##' Note that \code{\link{wl2i}} generates sequences of indices along the wavelength axis.
##' 
##' \code{seq} had to be implemented as S3 method as the generic has only \dots{} arguments (on which
##' no dispatch with differing types is possible).
##' 
##' \code{\link[base]{seq_along}} is not generic, but returns a sequence of the \code{length} of the
##' object. As hyperSpec provides a Method \code{\link{length}}, it can be used. The result is a
##' sequence of indices for the spectra.
##' 
##' @aliases seq seq,hyperSpec-method
##' @param x the hyperSpec object
##' @param from,to arguments handed to \code{\link[base]{seq.int}}
##' @param ... arguments for \code{\link[base]{seq}}, namely \code{by}, \code{length.out}
##' @param index should a vector with indices be returned rather than a hyperSpec object?
##' @return a numeric or hyperSpec object, depending on \code{index}.
##' @author C. Beleites
##' @seealso \code{\link{wl2i}} to construct sequences of wavelength indices.
##' 
##' \code{\link[base]{seq}}
##' @rdname seq
##' @export
##' @method seq hyperSpec
##' @keywords manip
##' @examples
##' 
##' seq (flu, index = TRUE)
##' seq_along (flu)
##' seq (flu, length.out = 3, index = TRUE) # return value is numeric, not integer! 
##' seq (flu, by = 2, index = TRUE) 		    # return value is numeric, not integer! 
##' 
##' plot (flu, col = "darkgray") 
##' plot (seq (flu, by = 2), add = TRUE, col= "red")
##' plot (seq (flu, length.out = 2), add = TRUE, col= "blue")
##' 
### needs to be an S3 function as S4 ... dispatch has to have the same signature for all parameters
seq.hyperSpec <- function (x, from = 1, to = nrow (x),
                           ..., index = FALSE){
  validObject (x)

  s <- seq (from = from, to = to, ...)

  if (index)
    s
  else {
    .extract (x, i = s)
	}
}

## internal abbreviation function

row.seq <- function (x, from = 1, to = nrow (x@data), ...){
  if (nrow (x@data) == 0)
    integer (0)
  else
    seq (from = from, to = to, ...)
}

col.seq <- function (x, from = 1, to = ncol (x@data), ...){
  if (ncol (x@data) == 0)
    integer (0)
  else
	seq (from = from, to = to, ...)
}

wl.seq <- function (x, from = 1, to = ncol (x@data$spc), ...){
  if (ncol (x@data$spc) == 0)
    integer (0)
  else
    seq (from = from, to = to, ...)
}
