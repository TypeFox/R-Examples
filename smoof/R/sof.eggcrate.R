#' Egg Crate Function
#'
#' This single-objective function follows the definition
#' \deqn{f(\mathbf{x}) = \mathbf{x}_1^2 + \mathbf{x}_2^2 + 25(\sin^2(\mathbf{x}_1) + \sin^2(\mathbf{x}_2))}
#' with \eqn{\mathbf{x}_i \in [-5, 5]} for \eqn{i = 1, 2}.
#'
#' @template ret_smoof_single
#' @export
makeEggCrateFunction = function() {
  makeSingleObjectiveFunction(
    name = "Egg Crate Function",
    id = "eggCrate_2d",
    fn = function(x) {
      assertNumeric(x, len = 2L, any.missing = FALSE, all.missing = FALSE)
      x[1]^2 + x[2]^2 + 25 * (sin(x[1])^2 + sin(x[2])^2)
    },
    par.set = makeNumericParamSet(
      len = 2L,
      id = "x",
      lower = c(-5, -5),
      upper = c(5, 5),
      vector = TRUE
    ),
    tags = attr(makeEggCrateFunction, "tags"),
    global.opt.params = c(0, 0),
    global.opt.value = 0
  )
}

class(makeEggCrateFunction) = c("function", "smoof_generator")
attr(makeEggCrateFunction, "name") = c("Egg Crate")
attr(makeEggCrateFunction, "type") = c("single-objective")
attr(makeEggCrateFunction, "tags") = c("single-objective", "continuous", "separable", "non-scalable")
