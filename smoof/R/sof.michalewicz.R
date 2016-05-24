#' Michalewicz Function
#'
#' Highly multimodal single-objective test function with \eqn{n!} local minima
#' with the formula:
#' \deqn{f(\mathbf{x}) = -\sum_{i=1}^{n} \sin(\mathbf{x}_i) \cdot \left(\sin\left(\frac{i \cdot \mathbf{x}_i}{\pi}\right)\right)^{2m}.}
#' The recommended value \eqn{m = 10}, which is used as a default in the
#' implementation.
#'
#' @note The location of the global optimum s varying based on both
#' the dimension and \eqn{m} parameter and is thus not provided in the
#' implementation.
#'
#' @references Michalewicz, Z.: Genetic Algorithms + Data Structures = Evolution
#' Programs. Berlin, Heidelberg, New York: Springer-Verlag, 1992.
#'
#' @template arg_dimensions
#' @param m [\code{integer(1)}]\cr
#'   \dQuote{Steepness} parameter.
#' @template ret_smoof_single
#' @export
makeMichalewiczFunction = function(dimensions, m = 10) {
  assertCount(dimensions, na.ok = FALSE)
  assertNumber(m, na.ok = FALSE)
  force(m)
  makeSingleObjectiveFunction(
    name = paste(dimensions, "-d Michalewicz Function (m = ", m, ")", sep = ""),
    id = paste0("michalewicz_", dimensions, "d_m", m),
    fn = function(x) {
      assertNumeric(x, len = 2L, any.missing = FALSE, all.missing = FALSE)
      i = 1:length(x)
      (-1) * sum(sin(x) * (sin((i * x^2) / pi)^(2 * m)))
    },
    par.set = makeNumericParamSet(
      len = 2L,
      id = "x",
      lower = c(0, 0),
      upper = c(pi, pi),
      vector = TRUE
    ),
    tags = attr(makeMichalewiczFunction, "tags"),
    global.opt.params = getMichalewiczGlobalOptimum(dimensions, m)$global.opt.params,
    global.opt.value = getMichalewiczGlobalOptimum(dimensions, m)$global.opt.value
  )
}

class(makeMichalewiczFunction) = c("function", "smoof_generator")
attr(makeMichalewiczFunction, "name") = c("Michalewicz")
attr(makeMichalewiczFunction, "type") = c("single-objective")
attr(makeMichalewiczFunction, "tags") = c("single-objective", "continuous", "multimodal", "scalable")

# Helper to determine global optimum based on dimension and m parameter.
getMichalewiczGlobalOptimum = function(dimensions = 2L, m = 10) {
  #FIXME: we know optimum only for m = 10
  if (m != 10) {
    return(NULL)
  }
  if (dimensions == 2L) {
    return(list(global.opt.params = c(2.20290552014618, 1.57079632677565), global.opt.value = -1.80130341009855321))
  }
  return(list())
}
