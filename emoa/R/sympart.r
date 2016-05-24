##
## sympart.r - SYM-PART test function from the CEC 2007 competition
##
## Author:
##  Olaf Mersmann (OME) <olafm@statistik.tu-dortmund.de>
##

##' Functions from the CEC 2007 EMOA competition.
##'
##' @param x Parmater vector.
##' @return Function value.
##'
##' @author Olaf Mersmann \email{olafm@@statistik.tu-dortmund.de}
##' @keywords optimize
##' @export
##' @rdname cec2007_functions
sympart <- function(x) {
  #stopifnot(length(x) >= 2)
  #stopifnot(length(x) %% 2 == 0)
  .Call(do_sympart, as.numeric(x))
}
