#' Generator for function with multiple peaks following the multiple peaks model 2.
#'
#' @param n.peaks [\code{integer(1)}]\cr
#'   Desired number of peaks, i. e., number of (local) optima.
#' @template arg_dimensions
#' @param topology [\code{character(1)}]\cr
#'   Type of topology. Possible values are \dQuote{random} and \dQuote{funnel}.
#' @param seed [\code{integer(1)}]\cr
#'   Seed for the random numbers generator.
#' @param rotated [\code{logical(1)}]\cr
#'   Should the peak shapes be rotated? This parameter is only relevant in case
#'   of elliptically shaped peaks.
#' @param peak.shape [\code{character(1)}]\cr
#'   Shape of peak(s). Possible values are \dQuote{ellipse} and \dQuote{sphere}.
#' @return [\code{smoof_single_objective_function}]
#' @examples
#' \dontrun{
#' fn = makeMPM2Function(n.peaks = 10L, dimensions = 2L,
#'   topology = "funnel", seed = 123, rotated = TRUE, peak.shape = "ellipse")
#' if (require(plot3D)) {
#'   plot3D(fn)
#' }
#' }
#' \dontrun{
#' fn = makeMPM2Function(n.peaks = 5L, dimensions = 2L,
#'   topology = "random", seed = 134, rotated = FALSE)
#' plot(fn, render.levels = TRUE)
#' }
#'
#' @references See the \href{https://ls11-www.cs.uni-dortmund.de/_media/techreports/tr15-01.pdf}{technical report}
#' of multiple peaks model 2 for an in-depth description of the underlying algorithm.
#'
#' @author \R interface by Jakob Bossek. Original python code provided by the Simon Wessing.
#'
#' @export
makeMPM2Function = function(n.peaks, dimensions, topology, seed, rotated = TRUE, peak.shape = "ellipse") {
  if (isWindows()) {
    stopf("No support for the multiple peaks model 2 generator at the moment.")
  }

  # do some sanity checks
  n.peaks = convertInteger(n.peaks)
  dimensions = convertInteger(dimensions)
  seed = convertInteger(seed)
  assertInt(n.peaks, lower = 1L, na.ok = FALSE)
  assertInt(dimensions, lower = 1L, na.ok = FALSE)
  assertChoice(topology, choices = c("random", "funnel"))
  assertInt(seed, lower = 1L, na.ok = FALSE)
  assertLogical(rotated, any.missing = FALSE)
  assertChoice(peak.shape, choices = c("ellipse", "sphere"))

  # touch vars
  force(n.peaks)
  force(dimensions)
  force(topology)
  force(seed)
  force(rotated)
  force(peak.shape)

  # build parameter set (bounds are [0, 1]^d)
  par.set = makeNumericParamSet("x", len = dimensions, lower = 0, upper = 1)

  # import rPython namespace
  BBmisc::requirePackages("_rPython", why = "smoof::makeMultiplePeaksModel2Function")

  # load funnel generator to global environemt
  eval(rPython::python.load(system.file("mpm2.py", package = "smoof")), envir = .GlobalEnv)

  local.opt.params = eval(rPython::python.call("getLocalOptimaParams", n.peaks, dimensions, topology, seed, rotated, peak.shape))
  if (n.peaks == 1)
    local.opt.params = list(local.opt.params)
  global.opt.params = eval(rPython::python.call("getGlobalOptimaParams", n.peaks, dimensions, topology, seed, rotated, peak.shape))

  smoof.fn = makeSingleObjectiveFunction(
    name = sprintf("Funnel_%i_%i_%i_%s_%s%s", n.peaks, dimensions, seed, topology, peak.shape, ifelse(rotated, "_rotated", "")),
    description = sprintf("Funnel-like function\n(n.peaks: %i, dimension: %i, topology: %s, seed: %i, rotated: %s, shape: %s)",
      n.peaks, dimensions, topology, seed, rotated, peak.shape),
    fn = function(x) {
      assertNumeric(x, len = dimensions, any.missing = FALSE, all.missing = FALSE)
      rPython::python.call("evaluateProblem", as.numeric(x), n.peaks, dimensions, topology, seed, rotated, peak.shape)
    },
    par.set = par.set,
    tags = c("non-separable", "scalable", "continuous", "multimodal")
  )
  smoof.fn = setAttribute(smoof.fn, "local.opt.params", local.opt.params)
  smoof.fn = setAttribute(smoof.fn, "local.opt.values", lapply(local.opt.params, smoof.fn))
  smoof.fn = setAttribute(smoof.fn, "global.opt.params", global.opt.params)
  smoof.fn = setAttribute(smoof.fn, "global.opt.value", smoof.fn(global.opt.params))
  return(smoof.fn)
}

class(makeMPM2Function) = c("function", "smoof_generator")
attr(makeMPM2Function, "name") = c("Multiple peals model 2 function generator")
attr(makeMPM2Function, "type") = c("single-objective")
