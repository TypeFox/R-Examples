##' colSums, colMeans, rowSums and rowMeans functions for hyperSpec objects
##'
##' hyperSpec objects can use the base functions \code{\link[base]{colMeans}},
##' \code{\link[base]{colSums}}, \code{\link[base]{rowMeans}} and \code{\link[base]{rowSums}}.
##'
##' @param x hyperSpec object
##' @param label.spc labels for the intensity axis for loadings-like (col) statistics
##' @param label.wavelength labels for the wavelength axis for scores-like (row) statistics
##' @param na.rm,... further parameters to the base functions
##'
##' \code{na.rm} defaults to \code{TRUE} for hyperSpec objects.
##' @seealso \link[base]{colSums}
##' @rdname colSums
##' @name colSums
NULL
 
##' @noRd
setGeneric ('colMeans')#, package = 'matrixStats')

##' @rdname colSums
##' @export
##' @examples
##' colMeans (flu)
 setMethod ("colMeans", signature = signature (x = "hyperSpec"), function (x, na.rm = TRUE, ..., label.spc){
   result <- colMeans (x@data$spc, na.rm = na.rm, ...)
   if (is.matrix (result) && ncol (result) != nwl (x) && nrow (result) == nwl (x))
      result <- t (result)

   decomposition (x, result, scores = FALSE, label.spc = label.spc)
}) 

.test (colMeans) <- function (){
  checkEqualsNumeric (colMeans (flu)[[]], colMeans (flu [[]]))
  
  checkEqualsNumeric (colMeans (fluNA)[[]], colMeans (fluNA [[]], na.rm = TRUE))
  checkEqualsNumeric (colMeans (fluNA, na.rm = FALSE)[[]], colMeans (fluNA [[]], na.rm = FALSE))
}

##' @noRd
setGeneric ('colSums') #, package = 'matrixStats')

##' @rdname colSums
##' @export
##' @examples
##' colSums (flu)
setMethod ("colSums", signature = signature (x = "hyperSpec"), function (x, na.rm = TRUE, ..., label.spc){
   result <- colSums (x@data$spc, na.rm = na.rm, ...)
   if (is.matrix (result) && ncol (result) != nwl (x) && nrow (result) == nwl (x))
      result <- t (result)

   decomposition (x, result, scores = FALSE, label.spc = label.spc)
}) 

.test (colSums) <- function (){
  checkEqualsNumeric (colSums (flu)[[]], colSums (flu [[]]))

  checkEqualsNumeric (colSums (fluNA)[[]], colSums (fluNA [[]], na.rm = TRUE))
  checkEqualsNumeric (colSums (fluNA, na.rm = FALSE)[[]], colSums (fluNA [[]], na.rm = FALSE))
}

##' @noRd
setGeneric ('rowMeans') #, package = 'matrixStats')

##' @rdname colSums
##' @export
##' @examples
##' colSums (flu)
setMethod ("rowMeans", signature = signature (x = "hyperSpec"), function (x, na.rm = TRUE, ..., label.wavelength){
   result <- rowMeans (x@data$spc, na.rm = na.rm, ...)
   if (is.matrix (result) && nrow (result) != nrow (x) && ncol (result) == nrow (x))
      result <- t (result)

   decomposition (x, result, scores = TRUE, label.wavelength = label.wavelength)
}) 

.test (rowMeans) <- function (){
  checkEqualsNumeric (rowMeans (flu)[[]], rowMeans (flu [[]]))

  checkEqualsNumeric (rowMeans (fluNA)[[]], rowMeans (fluNA [[]], na.rm = TRUE))
  checkEqualsNumeric (rowMeans (fluNA, na.rm = FALSE)[[]], rowMeans (fluNA [[]], na.rm = FALSE))
}

##' @noRd
setGeneric ('rowSums') #, package = 'matrixStats')

##' @rdname colSums
##' @export
##' @examples
##' rowSums (flu)
setMethod ("rowSums", signature = signature (x = "hyperSpec"), function (x, na.rm = TRUE, ..., label.wavelength){
   result <- rowSums (x@data$spc, na.rm = na.rm, ...)
   if (is.matrix (result) && nrow (result) != nrow (x) && ncol (result) == nrow (x))
      result <- t (result)

   decomposition (x, result, scores = TRUE, label.wavelength = label.wavelength)
}) 

.test (rowSums) <- function (){
  checkEqualsNumeric (rowSums (flu)[[]], rowSums (flu [[]]))

  checkEqualsNumeric (rowSums (fluNA)[[]], rowSums (fluNA [[]], na.rm = TRUE))
  checkEqualsNumeric (rowSums (fluNA, na.rm = FALSE)[[]], rowSums (fluNA [[]], na.rm = FALSE))
}

