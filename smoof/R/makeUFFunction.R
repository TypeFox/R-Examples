#' Generator for the functions UF1, ..., UF10 of the CEC 2009.
#'
#' @template arg_dimensions
#' @param id [\code{integer(1)}]\cr
#'   Instance identifier. Integer value between 1 and 10.
#' @return [\code{smoof_single_objective_function}]
#'
#' @author Jakob Bossek \email{j.bossek@@gmail.com}
#' @note The implementation is based on the original CPP implemenation by Qingfu
#' Zhang, Aimin Zhou, Shizheng Zhaoy, Ponnuthurai Nagaratnam Suganthany, Wudong Liu
#' and Santosh Tiwar.
#'
#' @export
makeUFFunction = function(dimensions, id) {
  # do some sanity checks
  dimensions = convertInteger(dimensions)
  id = convertInteger(id)
  #FIXME: I guess the lowest possible search space dimension is 3
  assertInt(dimensions, lower = 3L, na.ok = FALSE)
  assertInt(id, lower = 1L, upper = 10L, na.ok = FALSE)

  # touch vars
  force(dimensions)
  force(id)

  # problems 1, ..., 7 are 2D, the remainder is 3D
  n.objectives = 2L
  if (id > 7L) {
    n.objectives = 3L
  }

  makeMultiObjectiveFunction(
    name = sprintf("UF%i", id),
    id = paste0("UF", id, "_", dimensions, "d"),
    description = sprintf("One of the CEC 2009 functions."),
    fn = function(x) {
      .Call("evaluateUFFunction", as.integer(id), as.numeric(x), as.integer(dimensions))
    },
    par.set = makeUFParamSet(id, dimensions),
    n.objectives = n.objectives
  )
}

class(makeUFFunction) = c("function", "smoof_generator")
attr(makeUFFunction, "name") = c("UF1, ..., UF10 of the CEC 2009")
attr(makeUFFunction, "type") = c("Multi-objective")

# Helper function to determine the search space.
# Mostly [0, 1] x [-1, 1]^{n-1} for 2D (with two exceptions)
# and [0, 1]^2 x [-2, 2]^{n-2} for 3D
#
# @param id [integer(1)]
#   UF function id.
# @param dimension [integer(1)]
#   Problem space dimension.
# @return [ParamSet]
makeUFParamSet = function(id, dimensions) {
  n.objectives = if (id <= 7) 2L else 3L
  if (n.objectives == 2L) {
    if (id == 3L) {
      # UF3
      lower = rep(0, dimensions)
      upper = rep(1, dimensions)
    } else if (id == 4L) {
      # UF4
      lower = c(0, rep(-2, dimensions - 1L))
      upper = c(1, rep(2, dimensions - 1L))
    } else {
      # UF1, UF2, UF5, UF6, UF7
      lower = c(0, rep(-1, dimensions - 1L))
      upper = c(1, rep(1, dimensions - 1L))
    }
  } else {
    # UF8, UF9, UF10
    lower = c(0, 0, rep(-2, dimensions - 2L))
    upper = c(1, 1, rep(2, dimensions - 2L))
  }
  makeNumericParamSet("x", lower = lower, upper = upper, len = dimensions)
}
