#' Extracts ParamSet from soobench function.
#'
#' @param fn [\code{\link[soobench]{soo_function}}]\cr
#'   Source function.
#' @return [\code{\link{ParamSet}}]
#' @export

# FIXME: example disabled because of annying soobench interface change
# @examples
# library(soobench)
# fn = ackley_function(4)
# par.set = extractParamSetFromSooFunction(fn)
# print(par.set)
extractParamSetFromSooFunction = function(fn) {
  assertClass(fn, "soo_function")
  makeNumericParamSet(
    len = soobench::number_of_parameters(fn),
    id = "x",
    lower = soobench::lower_bounds(fn),
    upper = soobench::upper_bounds(fn),
    vector = FALSE
  )
}
