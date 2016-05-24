#' the function soft() encodes the soft-thresholding operator
#' @param x Float
#' @param d Float
#' @keywords internal
soft <-
function(x,d) return(sign(x)*pmax(0, abs(x)-d))
