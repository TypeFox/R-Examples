#' String representation of ecd
#' 
#' A string representation of an ecd object. Can be used for warning or error.
#'
#' @param object An object of ecd class
#' @param full logical, indicating if long form (multiple lines) should be rendered.
#'
#' @return character
#'
#' @keywords constructor
#'
#' @export
#'
#' @importFrom Rmpfr asNumeric
#'
#' @examples
#' ecd.toString(ecd(-1,1, sigma=0.1))
#'
### <======================================================================>
"ecd.toString" <- function(object, full=FALSE)
{
    a = sprintf("%.3f", object@alpha)
    r = sprintf("%.3f", object@gamma)
    s = sprintf("%e", asNumeric(object@sigma))
    m = sprintf("%e", asNumeric(object@mu))
    b = sprintf("%.2f", object@beta)
    l = sprintf("%.2f", object@lambda)
    R = sprintf("%.3f", object@R)
    t = sprintf("%.3f", object@theta)
    d = sprintf("%.1f", object@theta/pi*180)

    if (object@cusp > 0) {
        paste("ecd.cusp(",
              "a=", a, "r=", r, "s=", s,
              "m=", m, "b=", b,
              "cusp=", object@cusp, ")")
    } else if (object@lambda == 3) {
        paste("ecd(",
              "a=", a, "r=", r, "s=", s,
              "m=", m, "b=", b, "R=", R, "d=", d, ")")
    } else {
        paste("ecd.lambda(", "l=", l, "s=", s, "m=", m,"b=", b, ")")
    }
}
### <---------------------------------------------------------------------->

