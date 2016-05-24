#' Return ids of parameters in parameter set.
#'
#' Useful if vectors are included.
#' @template arg_parset
#' @param repeated [\code{logical(1)}]\cr
#'   Should ids be repeated length-times if parameter is a vector?
#'   Default is \code{FALSE}.
#' @param with.nr [\code{logical(1)}]\cr
#'   Should number from 1 to length be appended to id if \code{repeated} is \code{TRUE}?
#'   Otherwise ignored.
#'   Default is \code{FALSE}.
#' @return [\code{character}].
#' @export
#' @examples
#' ps = makeParamSet(
#'   makeNumericParam("u"),
#'   makeIntegerVectorParam("v", len = 2)
#' )
#' getParamIds(ps)
#' getParamIds(ps, repeated = TRUE)
#' getParamIds(ps, repeated = TRUE, with.nr = TRUE)
getParamIds = function(par.set, repeated = FALSE, with.nr = FALSE) {
  assertClass(par.set, "ParamSet")
  assertFlag(repeated)
  assertFlag(with.nr)
  getParamIds2(par.set$pars, repeated, with.nr)
}

getParamIds2 = function(pars, repeated = FALSE, with.nr = FALSE) {
  ns = lapply(pars, function(x) {
    if (repeated && x$type %in% c(
      "numericvector", "integervector", "discretevector",
      "logicalvector", "charactervector")) {
      n = x$len
      if (n > 1 && with.nr)
        paste(rep(x$id, n), 1:n, sep = "")
      else
        rep(x$id, n)
    } else {
      x$id
    }
  })
  as.character(do.call(c, ns))
}

