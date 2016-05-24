#' Generic Test
#' 
#' This generic function only exists to test that the rexygen2 parser work
#' correctly. Just ignore it.
#'
#' @param x Object
#' @param ... Object
#' @param methodParam Object
#' 
#' @include S4-generics.R
#' @rdname genericTest
#' @export
.genericTest(x, ...) %g% standardGeneric(".genericTest")

#' @export
#' @rdname genericTest
.genericTest(x ~ numeric, ..., methodParam = function() 1) %m% {
  x
}
