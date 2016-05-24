.sample <- function (x, size = nrow (x), replace = FALSE, prob = NULL) {
	validObject (x)
	
	s <- sample.int (nrow (x@data), size = size, replace = replace, prob = prob)
	
	x [s]
}

##' Random Samples and Permutations
##' Take a sample of the specified size from the elements of x with or without
##' replacement.
##' 
##' @rdname sample
##' @docType methods
##' @param x The hyperSpec object, data.frame or matrix to sample fromto sample from
##' @param size positive integer giving the number of spectra (rows) to choose.
##' @param replace Should sampling be with replacement?
##' @param prob A vector of probability weights for obtaining the elements of
##'   the vector being sampled.
##' @return a hyperSpec object, data.frame or matrix with \code{size} rows for \code{sample}, and an
##' integer vector for \code{isample} that is suitable for indexing (into the spectra) of x.
##' @author C. Beleites
##' @seealso \code{\link[base]{sample}}
##' @keywords methods distribution
##' @export
##' @examples
##' 
##' sample (flu, 3)
##' 
##' plot (flu, col = "darkgray")
##' plot (sample (flu, 3), col = "red", add = TRUE)
##' 
##' plot (flu, col = "darkgray")
##' plot (sample (flu, 3, replace = TRUE), col = "#0000FF80", add = TRUE,
##'       lines.args = list (lwd = 2));
##' 
setMethod ("sample", signature = signature (x = "hyperSpec"), .sample)

##' \code{isample} returns an vector of indices, \code{sample} returns the
##' corresponding hyperSpec object.
##'
##' @rdname sample
##' @return vector with indices suitable for row-indexing x
##' @export
##' @examples
##' isample (flu, 3)
##' isample (flu, 3, replace = TRUE)
##' isample (flu, 8, replace = TRUE)

isample <- function (x, size = nrow (x), replace = FALSE, prob = NULL) {
  chk.hy (x)
  validObject (x)

  sample.int (nrow (x), size = size, replace = replace, prob = prob)
}


##' @rdname sample
##' @param drop see \code{\link[base]{drop}}: by default, do not drop dimensions of the result
##' @export
##' @examples
##' sample (cars, 2)
setMethod ("sample", signature = signature (x = "data.frame"),
           function (x, size = nrow (x), replace = FALSE, prob = NULL, drop = FALSE) {
             x [sample.int (nrow (x), size = size, replace = replace, prob = prob), , drop = drop]
           }
           )
##' @rdname sample
##' @export
##' @examples
##' sample (matrix (1:24, 6), 2)
setMethod ("sample", signature = signature (x = "matrix"),
           function (x, size = nrow (x), replace = FALSE, prob = NULL, drop = FALSE) {
             x [sample.int (nrow (x), size = size, replace = replace, prob = prob), , drop = drop]
           }
           )
