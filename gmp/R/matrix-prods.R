if(FALSE) {## a nice idea --- but not working:  setOldClass() is fine,
    ## but then dispatch here only works if you create  new("bigz",...) objects
    ## (though dispatch *does* work for asNumeric(.)  -- really just R bug ???
    ## see also  'if(FALSE)'  in ./bigq.R

 ## NOTE: %*% is an S4, but *not* an S3 generic ==> Let's use S4 methods here

 setMethod("%*%", signature(x = "bigz", y = "bigz"),
           function(x,y) .Call(matrix_mul_z, x, y, FALSE, FALSE))
 setMethod("%*%", signature(x = "bigz", y = "ANY"),
           function(x,y) .Call(matrix_mul_z, x, y, FALSE, FALSE))
 setMethod("%*%", signature(x = "ANY", y = "bigz"),
           function(x,y) .Call(matrix_mul_z, x, y, FALSE, FALSE))
 setMethod("crossprod", signature(x = "bigz", y = "bigz"),
           function(x,y) .NotYetImplemented())
 setMethod("crossprod", signature(x = "bigz", y = "ANY"),
           function(x,y) .NotYetImplemented())
 setMethod("crossprod", signature(x = "ANY", y = "bigz"),
           function(x,y) .NotYetImplemented())
 setMethod("tcrossprod", signature(x = "bigz", y = "bigz"),
           function(x,y) .NotYetImplemented())
 setMethod("tcrossprod", signature(x = "bigz", y = "ANY"),
           function(x,y) .NotYetImplemented())
 setMethod("tcrossprod", signature(x = "ANY", y = "bigz"),
           function(x,y) .NotYetImplemented())

 setMethod("%*%", signature(x = "bigq", y = "bigq"),
           function(x,y) .Call(matrix_mul_q, x, y, FALSE, FALSE))
 setMethod("%*%", signature(x = "bigq", y = "ANY"),
           function(x,y) .Call(matrix_mul_q, x, y, FALSE, FALSE))
 setMethod("%*%", signature(x = "ANY", y = "bigq"),
           function(x,y) .Call(matrix_mul_q, x, y, FALSE, FALSE))
 setMethod("crossprod", signature(x = "bigq", y = "bigq"),
           function(x,y) .NotYetImplemented())
 setMethod("crossprod", signature(x = "bigq", y = "ANY"),
           function(x,y) .NotYetImplemented())
 setMethod("crossprod", signature(x = "ANY", y = "bigq"),
           function(x,y) .NotYetImplemented())
 setMethod("tcrossprod", signature(x = "bigq", y = "bigq"),
           function(x,y) .NotYetImplemented())
 setMethod("tcrossprod", signature(x = "bigq", y = "ANY"),
           function(x,y) .NotYetImplemented())
 setMethod("tcrossprod", signature(x = "ANY", y = "bigq"),
           function(x,y) .NotYetImplemented())

} else { ## less nice -- S3-only -- way: --------------------------------------------
 `%*%` <- function(x,y) UseMethod("%*%")

 `%*%.default` <- function(x,y) {
     ## dispatch on y (!)
     if(inherits(y, "bigz"))
         .Call(matrix_mul_z, x, y, 0L)
     else if(inherits(y, "bigq"))
         .Call(matrix_mul_q, x, y, 0L)
     else base::"%*%"(x,y)
 }

 ## matrix_mul_z (SEXP a, SEXP b, SEXP right, SEXP trans)
 ##     right if(right), compute b %*% T(a) , else T(a) %*% b
 ##     trans if(trans), T(a) := t(a) = a',  else T(a) := a

 `%*%.bigz` <- function(x,y) .Call(matrix_mul_z, x, y, 0L)
 `%*%.bigq` <- function(x,y) .Call(matrix_mul_q, x, y, 0L)

 crossprod  <- function(x,y=NULL) UseMethod("crossprod")
 tcrossprod <- function(x,y=NULL) UseMethod("tcrossprod")

 crossprod.default <- function(x,y=NULL) {
     if(is.null(y))
	 return(base::crossprod(x))
     if(inherits(y, "bigz"))
	 .Call(matrix_mul_z, x, y, 1L)
     else if(inherits(y, "bigq"))
	 .Call(matrix_mul_q, x, y, 1L)
     else base::crossprod(x,y)
 }

 crossprod.bigz <- function(x,y=NULL) {
     if(is.null(y))
	 .Call(matrix_crossp_z, x, FALSE)
     else if(inherits(y, "bigq"))
	 .Call(matrix_mul_q, x, y, 1L)
     else
	 .Call(matrix_mul_z, x, y, 1L)
 }

 crossprod.bigq <- function(x,y=NULL)
 {
     if(is.null(y))
	 .Call(matrix_crossp_q, x, FALSE)
     else
	 .Call(matrix_mul_q, x, y, 1L)
 }

 ##-----------------------------------------------

 tcrossprod.default <- function(x,y=NULL) {
     if(is.null(y))
	 return(base::tcrossprod(x))
     if(inherits(y, "bigz"))
	 .Call(matrix_mul_z, x, y, 2L)
     else if(inherits(y, "bigq"))
	 .Call(matrix_mul_q, x, y, 2L)
     else base::tcrossprod(x,y)
 }

 tcrossprod.bigz <- function(x,y=NULL) {
     if(is.null(y))
	 .Call(matrix_crossp_z, x, TRUE)
     else if(inherits(y, "bigq"))
	 .Call(matrix_mul_q, x, y, 2L)
     else
	 .Call(matrix_mul_z, x, y, 2L)
 }

 tcrossprod.bigq <- function(x,y=NULL) {
     if(is.null(y))
	 .Call(matrix_crossp_q, x, TRUE)
     else
	 .Call(matrix_mul_q, x, y, 2L)
 }

}# end{S3-way}
