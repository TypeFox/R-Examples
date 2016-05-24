#' The ecldOrEcd class
#'
#' The S4 class union of ecld and ecd, primarily used to define slot in \code{ecop.opt} class.
#' Its usage is rather cumbersome, so the end user should avoid it as much as possible.
#' 
#' @name ecldOrEcd-class
#' @include ecld-class.R
#' @include ecd-class.R
#'
#' @exportClass ecldOrEcd
setClassUnion("ecldOrEcd", c("ecld", "ecd"))

# end

