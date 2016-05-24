##' Arithmetical Operators: +, -, *, /, ^, \%\%, \%/\%, \%*\% for hyperSpec objects
##'
##' The arithmetical operators \code{+}, \code{-}, \code{*}, \code{/}, \code{\^}, \code{\%\%},
##' \code{\%/\%},  and \code{\%*\%} for \code{hyperSpec} objects. 
##' 
##' You can use these operators in different ways:
##' \preformatted{
##' e1 + e2
##' `+` (e1, e2)
##' 
##' x \%*\% y
##' `\%*\%`(x, y)
##' 
##' -x }
##' The arithmetical operators \code{+}, \code{-}, \code{*}, \code{/}, \code{^}, \code{\%\%},
##' \code{\%/\%}, and \code{\%*\%} work on the  spectra matrix of the \code{hyperSpec} object. They
##' have their usual meaning (see \code{\link[base]{Arithmetic}}).  The operators work also with one
##' \code{hyperSpec} object and a numeric object or a matrices of the same size as the spectra matrix
##' of the \code{hyperSpec} object.
##' 
##' With numeric vectors \code{\link[hyperSpec]{sweep}} is most probably more appropriate.
##'   
##' If you want to calculate on the extra data as well, use the data.frame \code{hyperSpec@@data}
##' directly or \code{\link[hyperSpec]{as.data.frame} (x)}.
##' @author C. Beleites
##' @title Arithmetical Operators for hyperSpec objects
##' @name Arith
##' @rdname Arith
##' @docType methods
##' @aliases Arith Arith,hyperSpec-method Arith,hyperSpec,hyperSpec-method
##' +,hyperSpec,hyperSpec-method -,hyperSpec,hyperSpec-method *,hyperSpec,hyperSpec-method
##' ^,hyperSpec,hyperSpec-method %%,hyperSpec,hyperSpec-method %/%,hyperSpec,hyperSpec-method
##' /,hyperSpec,hyperSpec-method
##' Arith,hyperSpec,matrix-method
##' Arith,hyperSpec,numeric-method
##' Arith,hyperSpec,missing-method
##' Arith,matrix,hyperSpec-method
##' Arith,numeric,hyperSpec-method
##' %*% %*%,hyperSpec,hyperSpec-method %*%,matrix,hyperSpec-method %*%,hyperSpec,matrix-method
##' @param e1,e2 or
##' @param x,y either two \code{hyperSpec} objects or
##'
##' one \code{hyperSpec} object and  matrix of same size as \code{hyperSpec[[]]} or
##'
##' a vector which length equalling either the number of rows or the number of wavelengths
##' of the hyperSpec object, or
##' 
##' a scalar (numeric of length 1).
##' @return \code{hyperSpec} object with the new spectra matrix.
##' @export
##' @keywords methods arith
##' @include paste.row.R
##' @include hyperspec-package.R
##' @include hyperspec-class.R
##' @concept hyperSpec arithmetic
##' @concept hyperSpec arithmetical operators
##' @concept hyperSpec plus
##' @concept hyperSpec division
##' @concept hyperSpec spectra conversion
##' @seealso
##' \code{\link[hyperSpec]{sweep-methods}} for calculations involving a vector and
##'   the spectral matrix.
##'   
##'   \code{\link[methods]{S4groupGeneric}} for group generic methods.
##' 
##'   \code{\link[base]{Arithmetic}} for the base arithmetic functions.
##' 
##'   \code{\link[hyperSpec]{Comparison}} for comparison operators, 
##'   \code{\link[hyperSpec]{Math}} for mathematical group generic 
##'   functions (Math and Math2 groups) working on \code{hyperSpec} objects.
##' @examples
##' flu + flu
##' 1 / flu
##' all((flu + flu - 2 * flu)[[]] == 0)
##' -flu
##' flu / flu$c


setMethod ("Arith", signature (e1 = "hyperSpec", e2 = "hyperSpec"),
           function (e1, e2){
             validObject (e1)
             validObject (e2)

             e1 <- .expand (e1, dim (e2) [c (1, 3)])
             e2 <- .expand (e2, dim (e1) [c (1, 3)])
             
             e1 [[]] <- callGeneric (e1[[]], e2[[]])
             e1
           }
           )

.arithx <- function (e1, e2){
  validObject (e1)

  if (missing (e2)){
    e1  [[]] <- callGeneric (e1 [[]])
    e1
  } else {
    e2 <- as.matrix (e2)

    ## called /only/ with e1 hyperSpec but e2 matrix-like
    e1 <- .expand (e1, dim (e2))
    e2 <- .expand (e2, dim (e1) [c (1, 3)])
    
    e1  [[]] <- callGeneric (e1 [[]], e2)
    e1
  }
}
##' @rdname Arith
setMethod ("Arith", signature (e1 = "hyperSpec", e2 = "numeric"), .arithx)
##' @rdname Arith
setMethod ("Arith", signature (e1 = "hyperSpec", e2 = "matrix"), .arithx)
##' @rdname Arith
setMethod ("Arith", signature (e1 = "hyperSpec", e2 = "missing"), .arithx)

.arithy <- function (e1, e2){
  e1 <- as.matrix (e1)
  validObject (e2)
  
  ## called /only/ with e2 hyperSpec but e1 matrix-like
  e1 <- .expand (e1, dim (e2) [c (1, 3)])
  e2 <- .expand (e2, dim (e1))
    
  e2  [[]] <- callGeneric (e1, e2 [[]])
  e2
}
##' @rdname Arith
setMethod ("Arith", signature (e1 = "numeric", e2 = "hyperSpec"), .arithy)
##' @rdname Arith
setMethod ("Arith", signature (e1 = "matrix", e2 = "hyperSpec"), .arithy)

##' @param m matrix
##' @param target.dim target size to expand the vector to for the sweep-shortcuts
##' @noRd
.expand <- function (m, target.dim) {
  m.dim = dim (m)
   
  if (m.dim [1]  == 1L & target.dim [1] > 1L)
      m <- m [rep (1, target.dim [1]),, drop = FALSE]

  if (is (m, "hyperSpec")) {
    if (m.dim [3] == 1L & target.dim [2] > 1L)
        m <- m [,, rep (1, target.dim [2]), wl.index = TRUE]
  } else {
    if (m.dim [2] == 1L & target.dim [2] > 1L)
        m <- m [, rep (1, target.dim [2]), drop = FALSE]
  }
  
  m
}


##' @rdname Arith
##' @concept hyperSpec matrix multiplication
##' @aliases %*% %*%,hyperSpec,hyperSpec-method %*%,matrix,hyperSpec-method
##' %*%,hyperSpec,matrix-method
##' @export "%*%"
##' @seealso  \code{\link[base]{matmult}} for matrix multiplications with \code{\%*\%}.
setMethod ("%*%", signature (x = "hyperSpec", y = "hyperSpec"),
           function (x, y){
             validObject (x)
             validObject (y)

             if (ncol(y) > 1)
               warning(paste("Dropping column(s) of y:", paste(colnames(y$..),
                                                               collapse = ", ")))

             x@data$spc <-  x@data$spc %*% y@data$spc
             .wl (x) <- y@wavelength
             x@label$.wavelength = y@label$.wavelength

             x
           }
           )

##' @rdname Arith
setMethod ("%*%", signature (x = "hyperSpec", y = "matrix"),
           function (x, y){
             validObject (x)
             x@data$spc <-  x@data$spc %*% y
             .wl (x) <- seq_len (ncol (y))
             x@label$.wavelength = NA
             x
           }
           )

##' @rdname Arith
setMethod ("%*%", signature (x = "matrix", y = "hyperSpec"),
           function (x, y){
             validObject (y)

             if (ncol(y) > 1)
               warning(paste("Dropping column(s) of y:", paste(colnames(y$..),
                                                               collapse = ", ")))
             y <- new ("hyperSpec",
                       wavelength = y@wavelength,
                       spc = x %*% y@data$spc,
                       log = y@log
                       )

             y
           }
           )

.test (Arith) <- function (){
  checkEqualsNumeric (as.matrix (flu - flu), rep (0, nrow (flu) * nwl (flu)))
  checkEqualsNumeric (as.matrix (flu - flu [1]), as.matrix (sweep (flu, 2, flu [1], `-`)))
  checkEqualsNumeric (as.matrix (flu - flu [,, 450]), as.matrix (sweep (flu, 1, flu [,, 450], `-`)))

  checkEqualsNumeric (as.matrix (flu / flu), rep (1, nrow (flu) * nwl (flu)))
  checkEqualsNumeric (as.matrix (flu / flu [1]), as.matrix (sweep (flu, 2, flu [1], `/`)))
  checkEqualsNumeric (as.matrix (flu / flu [,, 450]), as.matrix (sweep (flu, 1, flu [,, 450], `/`)))

  checkEqualsNumeric (as.matrix (flu + 1), as.matrix (flu) + 1) 
  checkEqualsNumeric (as.matrix (1 + flu), as.matrix (flu) + 1) 

  checkEqualsNumeric (as.matrix (-flu), - as.matrix (flu))
}
