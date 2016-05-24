#' Generator for the noiseless function set of the real-parameter Black-Box
#' Optimization Benchmarking (BBOB).
#'
#' @note It is possible to pass a matrix of parameters to the functions, where
#' each column consists of one parameter setting.
#'
#' @param dimension [\code{integer(1)}]\cr
#'   Problem dimension. Integer value between 2 and 40.
#' @param fid [\code{integer(1)}]\cr
#'   Function identifier. Integer value between 1 and 24.
#' @param iid [\code{integer(1)}]\cr
#'   Instance identifier. Integer value greater than or equal 1.
#' @return [\code{smoof_single_objective_function}]
#' @examples
#' # get the first instance of the 2D Sphere function
#' fn = makeBBOBFunction(dimension = 2L, fid = 1L, iid = 1L)
#' if (require(plot3D)) {
#'   plot3D(fn, contour = TRUE)
#' }
#' @references See the \href{http://coco.gforge.inria.fr/doku.php?id=bbob-2009-downloads}{BBOB website}
#' for a detailed description of the BBOB functions.
#' @export
makeBBOBFunction = function(dimension, fid, iid) {
  # do some sanity checks
  dimension = convertInteger(dimension)
  fid = convertInteger(fid)
  iid = convertInteger(iid)
  assertInt(dimension, lower = 2L, upper = 40L, na.ok = FALSE)
  assertInt(fid, lower = 1L, upper = 24L, na.ok = FALSE)
  assertInt(iid, lower = 1L, na.ok = FALSE)

  # touch vars
  force(dimension)
  force(fid)
  force(iid)

  # build parameter set (bounds are [-5, 5] for all BBOB funs)
  par.set = makeNumericParamSet("x", len = dimension, lower = -5, upper = 5)

  # get optimal values
  optimals = getOptimumForBBOBFunction(dimension, fid, iid)

  # get metadata, i. e., tags and name
  meta = mapBBOBFidToMetaData(fid)

  makeSingleObjectiveFunction(
    name = sprintf("BBOB_%i_%i_%i", dimension, fid, iid),
    description = sprintf("%i-th noiseless BBOB function\n(FID: %i, IID: %i, DIMENSION: %i)",
      fid, fid, iid, dimension),
    fn = function(x) {
      .Call("evaluateBBOBFunctionCPP", dimension, fid, iid, x)
    },
    par.set = par.set,
    tags = meta$tags,
    # all BBOB functions are vectorized
    vectorized = TRUE,
    global.opt.params = as.numeric(optimals[[1L]]),
    global.opt.value = as.numeric(optimals[[2L]])
  )
}

class(makeBBOBFunction) = c("function", "smoof_generator")
attr(makeBBOBFunction, "name") = c("Set of noiseless BOBB Function(s)")
attr(makeBBOBFunction, "type") = c("single-objective")

mapBBOBFidToMetaData = function(fid) {
  mapping = list(
    "1" = list(name = "Sphere", tags = c("unimodal", "separable", "differentiable", "continuous", "convex")),
    "2" = list(name = "Ellipsoidal", tags = c("unimodal", "separable", "differentiable", "continuous", "convex")),
    "3" = list(name = "Rastrigin", tags = c("multimodal", "separable", "differentiable", "continuous")),
    "4" = list(name = "Bueche-Rastrigin", tags = c("multimodal", "separable", "differentiable", "continuous")),
    "5" = list(name = "Linear Slope", tags = c("unimodal", "separable", "differentiable", "continuous")),
    "6" = list(name = "Attractive Sector", tags = c("unimodal", "continuous", "moderate-conditioned", "non-separable")),
    "7" = list(name = "Step Ellipsoidal", tags = c("unimodal", "moderate-conditioned", "non-separable", "non-differentiable")),
    "8" = list(name = "Rosenbrock (original)", tags = c("continuous", "differentiable", "non-separable", "scalable", "multimodal")),
    "9" = list(name = "Rosenbrock (rotated)", tags = c("continuous", "differentiable", "non-separable", "scalable", "multimodal")),
    "10" = list(name = "Ellipsoidal", tags = c("continuous", "differentiable", "non-separable", "scalable", "unimodal", "highly-conditioned")),
    "11" = list(name = "Discus", tags = c("continuous", "differentiable", "non-separable", "scalable", "unimodal", "highly-conditioned")),
    "12" = list(name = "Bent Cigar", tags = c("continuous", "differentiable", "non-separable", "scalable", "unimodal", "highly-conditioned")),
    "13" = list(name = "Sharp Ridge", tags = c("continuous", "differentiable", "non-separable", "scalable", "unimodal", "highly-conditioned")),
    "14" = list(name = "Different Powers", tags = c("continuous", "differentiable", "non-separable", "scalable", "unimodal", "highly-conditioned")),
    "15" = list(name = "Rastrigin", tags = c("continuous", "differentiable", "non-separable", "scalable", "multimodal", "adequate-global-structure")),
    "16" = list(name = "Weierstrass", tags = c("continuous", "differentiable", "non-separable", "scalable", "multimodal", "adequate-global-structure")),
    "17" = list(name = "Schaffers F7", tags = c("continuous", "differentiable", "non-separable", "scalable", "multimodal", "adequate-global-structure")),
    "18" = list(name = "Schaffers F7 (moderately ill-conditioned)", c("continuous", "differentiable", "non-separable", "scalable", "multimodal", "adequate-global-structure", "moderate-conditioned")),
    "19" = list(name = "Composite Griewank-Rosenbrock F8F2", c("continuous", "differentiable", "non-separable", "scalable", "multimodal", "adequate-global-structure")),
    "20" = list(name = "Schwefel", c("continuous", "differentiable", "non-separable", "scalable", "multimodal", "weak-global-structure")),
    "21" = list(name = "Gallagher's Gaussian 101-me Peaks", c("continuous", "differentiable", "non-separable", "scalable", "multimodal", "weak-global-structure")),
    "22" = list(name = "Gallagher's Gaussian 21-hi Peaks", c("continuous", "differentiable", "non-separable", "scalable", "multimodal", "weak-global-structure")),
    "23" = list(name = "Katsuura", c("continuous", "differentiable", "non-separable", "scalable", "multimodal", "weak-global-structure")),
    "24" = list(name = "Lunacek bi-Rastrigin", c("continuous", "differentiable", "non-separable", "scalable", "multimodal", "weak-global-structure"))
  )
  mapping[[fid]]
}

# Get the optimal parameter values and the optimal function value for a BBOB
# function.
getOptimumForBBOBFunction = function(dimension, fid, iid) {
  .Call("getOptimumForBBOBFunctionCPP", dimension, fid, iid)
}
