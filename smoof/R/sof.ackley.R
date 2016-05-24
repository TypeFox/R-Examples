#' Ackley Function
#'
#' Also known as \dQuote{Ackley's Path Function}.
#' Multimodal test function with its global optimum in the center of the defintion
#' space. The implementation is based on the formula
#' \deqn{f(\mathbf{x}) = -a \cdot \exp\left(-b \cdot \sqrt{\left(\frac{1}{n} \sum_{i=1}^{n} \mathbf{x}_i\right)}\right) - \exp\left(\frac{1}{n} \sum_{i=1}^{n} \cos(c \cdot \mathbf{x}_i)\right),}
#' with \eqn{a = 20}, \eqn{b = 0.2} and \eqn{c = 2\pi}. The feasible region is
#' given by the box constraints \eqn{\mathbf{x}_i \in [-32.768, 32.768]}.
#'
#' @references Ackley, D. H.: A connectionist machine for genetic hillclimbing.
#' Boston: Kluwer Academic Publishers, 1987.
#'
#' @template arg_dimensions
#' @template ret_smoof_single
#' @export
makeAckleyFunction = function(dimensions) {
  assertCount(dimensions, na.ok = FALSE)
  force(dimensions)
  makeSingleObjectiveFunction(
    name = paste(dimensions, "-d Ackley Function", sep = ""),
    id = paste0("ackley_", dimensions, "d"),
    fn = function(x) {
      assertNumeric(x, len = dimensions, any.missing = FALSE, all.missing = FALSE)
      n = length(x)
      a = 20
      b = 0.2
      c = 2 * pi
      d = mean(x^2)
      e = mean(cos(c * x))
      -a * exp(-b * sqrt(d)) - exp(e) + a + exp(1)
    },
    par.set = makeNumericParamSet(
      len = dimensions,
      id = "x",
      lower = -32.768,
      upper = 32.768,
      vector = TRUE
    ),
    tags = attr(makeAckleyFunction, "tags"),
    global.opt.params = rep(0, dimensions),
    global.opt.value = 0L
  )
}

class(makeAckleyFunction) = c("function", "smoof_generator")
attr(makeAckleyFunction, "name") = c("Ackley")
attr(makeAckleyFunction, "type") = c("single-objective")
attr(makeAckleyFunction, "tags") = c("single-objective", "continuous", "multimodal", "differentiable", "non-separable", "scalable")
