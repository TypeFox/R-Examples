##' deconstruct call to modifier into list
##' 
##' deconstructs the call to a modifier into a list for further
##' processing with the list-based parameter selection functions
##' @keywords internal
make_range_list <- function(call, method) {
  par_months <- eval(call$.months)
  par_variables <- eval(call$.variables)
  if (is.null(par_months)) {
    par_months <- -6:9
  }
  if (is.null(par_variables)) {
    selection <- list(method, par_months)
  } else {
    selection <- list(method, par_months, par_variables)
  }
  class(selection) <- c("tc_paramlist", "list")
  selection
}

##' Modifiers for climate parameter selection
##' 
##' These modifiers are used to select specific months from specific
##' climate parameters, and potentially transform the selections into
##' their respective sums or means. The modifiers can be chained
##' together using '+'.
##' @rdname treeclim-modifiers
##' @param .months numeric identifiers for the months (-1 for previous
##' January until 12 for current December, with -6 for previous June,
##' etc.)
##' @param .variables names of the variables the modifier shall be applied to
##' @export
.range <- function(.months = NULL, .variables = NULL) {
  make_range_list(match.call(), "full")
}

##' @rdname treeclim-modifiers
##' @export
.mean <- function(.months = NULL, .variables = NULL) {
  make_range_list(match.call(), "mean")
}

##' @rdname treeclim-modifiers
##' @export
.sum <- function(.months = NULL, .variables = NULL) {
  make_range_list(match.call(), "sum")
}
