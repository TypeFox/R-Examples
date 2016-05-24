##
## cec2009.r - Test functions from the CEC2009 competition
##
## Author:
##  Olaf Mersmann (OME) <olafm@statistik.tu-dortmund.de>
##

##' Functions from the CEC 2009 EMOA competition.
##'
##' @param x Parmater vector.
##' @return Function value.
##'
##' @author Olaf Mersmann \email{olafm@@statistik.tu-dortmund.de}
##' @keywords optimize
##' @export
##' @rdname cec2009_functions
UF1 <- function(x) {
  stopifnot(length(x) >= 3)
  .Call("do_UF1", as.numeric(x), PACKAGE="emoa")
}

##' @export
##' @rdname cec2009_functions
UF2 <- function(x) {
  stopifnot(length(x) >= 3)
  .Call("do_UF2", as.numeric(x), PACKAGE="emoa")
}

##' @export
##' @rdname cec2009_functions
UF3 <- function(x) {
  stopifnot(length(x) >= 3)
  .Call("do_UF3", as.numeric(x), PACKAGE="emoa")
}

##' @export
##' @rdname cec2009_functions
UF4 <- function(x) {
  stopifnot(length(x) >= 3)
  .Call("do_UF4", as.numeric(x), PACKAGE="emoa")
}

##' @export
##' @rdname cec2009_functions
UF5 <- function(x) {
  stopifnot(length(x) >= 3)
  .Call("do_UF5", as.numeric(x), PACKAGE="emoa")
}

##' @export
##' @rdname cec2009_functions
UF6 <- function(x) {
  stopifnot(length(x) >= 3)
  .Call("do_UF6", as.numeric(x), PACKAGE="emoa")
}

##' @export
##' @rdname cec2009_functions
UF7 <- function(x) {
  stopifnot(length(x) >= 3)
  .Call("do_UF7", as.numeric(x), PACKAGE="emoa")
}

##' @export
##' @rdname cec2009_functions
UF8 <- function(x) {
  stopifnot(length(x) >= 5)
  .Call("do_UF8", as.numeric(x), PACKAGE="emoa")
}

##' @export
##' @rdname cec2009_functions
UF9 <- function(x) {
  stopifnot(length(x) >= 5)
  .Call("do_UF9", as.numeric(x), PACKAGE="emoa")
}

##' @export
##' @rdname cec2009_functions
UF10 <- function(x) {
  stopifnot(length(x) >= 5)
  .Call("do_UF10", as.numeric(x), PACKAGE="emoa")
}
