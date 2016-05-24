#' Sample a random value from a parameter or a parameter set uniformly.
#'
#' @template desc_dep_pars_na
#'
#' @template arg_par_or_set
#' @template arg_disc_names
#' @template arg_trafo
#' @return The return type is determined by the type of the parameter. For a set a named list
#'   of such values in the correct order is returned.
#' @export
#' @examples
#' # bounds are necessary here, can't sample with Inf bounds:
#' u = makeNumericParam("x", lower = 0, upper = 1)
#' # returns a random number between 0 and 1:
#' sampleValue(u)
#'
#' p = makeDiscreteParam("x", values = c("a","b","c"))
#' # can be either "a", "b" or "c"
#' sampleValue(p)
#'
#' p = makeIntegerVectorParam("x", len = 2, lower = 1, upper = 5)
#' # vector of two random integers between 1 and 5:
#' sampleValue(p)
#'
#' ps = makeParamSet(
#'   makeNumericParam("x", lower = 1, upper = 10),
#'   makeIntegerParam("y", lower = 1, upper = 10),
#'   makeDiscreteParam("z", values = 1:2)
#' )
#' sampleValue(ps)
sampleValue = function(par, discrete.names = FALSE, trafo = FALSE) {
  UseMethod("sampleValue")
}

#' @export
sampleValue.Param = function(par, discrete.names = FALSE, trafo = FALSE) {
  type = par$type
  if (par$type %in% c("numeric", "numericvector", "integer", "integervector"))
    if (any(is.infinite(c(par$lower, par$upper))))
      stop("Cannot sample with Inf bounds!")
  if (!is.null(par$len) && is.na(par$len))
    stop("Cannot sample with NA length!")
  if (type == "numeric") {
    x = runif(1, min = par$lower, max = par$upper)
  } else if (type == "numericvector") {
    x = runif(par$len, min = par$lower, max = par$upper)
  } else if (type == "integer") {
    x = as.integer(round(runif(1, min = par$lower-0.5, max = par$upper+0.5)))
  } else if (type == "integervector") {
    x = as.integer(round(runif(par$len, min = par$lower-0.5, max = par$upper+0.5)))
  } else if (type %in% c("logical", "logicalvector")) {
    x = sample(c(TRUE, FALSE), par$len, replace = TRUE)
  } else if (type %in% c("discrete", "discretevector")) {
    x = sample(names(par$values), par$len, replace = TRUE)
    if (!discrete.names) {
      x = if (type  == "discretevector")
        par$values[x]
      else
        par$values[[x]]
    }
  } else if (type %in% c("function", "untyped", "character", "charactervector")) {
    stopf("Cannot generate random value for %s variable!", type)
  }
  if (trafo && !is.null(par$trafo))
    x = par$trafo(x)
  # if the components have names, set them
  if (!is.null(par$cnames))
    names(x) = par$cnames
  return(x)
}

#' @export
sampleValue.ParamSet = function(par, discrete.names = FALSE, trafo = FALSE) {
  # sample value for each param, do it until we a get one which is not forbidden
  repeat {
    val = lapply(par$pars, sampleValue, discrete.names = discrete.names, trafo = trafo)
    if (is.null(par$forbidden) || !isForbidden(par, val))
      break
  }
  # set conditional params to NA is condition not OK
  val = lapply(seq_along(val), function(i) {
    if (!is.null(par$pars[[i]]$requires) && !requiresOk(par, val, i)) {
      type = par$pars[[i]]$type
      type = switch(type,
        numericvector = "numeric",
        integervector = "integer",
        logicalvector = "logical",
        discrete = "character",
        discretevector = "character",
        type
      )
      as(NA, type)
     } else {
      val[[i]]
     }
  })
  names(val) = names(par$pars)
  return(val)
}


#' Sample n random values from a parameter or a parameter set uniformly.
#'
#' @template desc_dep_pars_na
#'
#' @template arg_par_or_set
#' @param n [\code{integer(1)}]\cr
#'   Number of values.
#' @template arg_disc_names
#' @template arg_trafo
#' @return [\code{list}]. For consistency always a list is returned.
#' @export
#' @examples
#' p = makeIntegerParam("x", lower = -10, upper = 10)
#' sampleValues(p, 4)
#'
#' p = makeNumericParam("x", lower = -10, upper = 10)
#' sampleValues(p, 4)
#'
#' p = makeLogicalParam("x")
#' sampleValues(p, 4)
#'
#' ps = makeParamSet(
#'   makeNumericParam("u", lower = 1, upper = 10),
#'   makeIntegerParam("v", lower = 1, upper = 10),
#'   makeDiscreteParam("w", values = 1:2)
#' )
#' sampleValues(ps, 2)
sampleValues = function(par, n, discrete.names = FALSE, trafo = FALSE) {
  assert(checkClass(par, "Param"), checkClass(par, "ParamSet"))
  n = asInt(n)
  assertFlag(discrete.names)
  replicate(n, sampleValue(par, discrete.names = discrete.names, trafo = trafo), simplify = FALSE)
}

