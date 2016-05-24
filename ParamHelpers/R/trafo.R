#' Transform a value.
#'
#' Transform a value with associated transformation function(s).
#' @template arg_par_or_set
#' @param x [any] \cr
#'   Single value to check.
#'   For a parameter set this must be a list in the correct order.
#' @return Transformed value.
#' @export
#' @examples
#' # transform simple parameter:
#' p = makeNumericParam(id="x", trafo=function(x) x^2)
#' trafoValue(p, 2)
#' # for a parameter set different transformation functions are possible:
#' ps = makeParamSet(
#'   makeIntegerParam("u", trafo=function(x) 2*x),
#'   makeNumericVectorParam("v", len=2, trafo=function(x) x/sum(x)),
#'   makeDiscreteParam("w", values=c("a", "b"))
#' )
#' # now the values of "u" and "v" are transformed:
#' trafoValue(ps, list(3, c(2, 4), "a"))
trafoValue = function(par, x) {
  if (inherits(par, "ParamSet"))
    Map(trafoValue, par$pars, x)
  else
    if(is.null(par$trafo))
      x
    else
      par$trafo(x)
}

#' Transform optimization path.
#'
#' Transform optimization path with associated transformation functions of parameters.
#' Can only be done when x values where added \dQuote{untransformed}.
#'
#' @param opt.path [\code{\link{OptPath}}]\cr
#'   Optimization path.
#' @return [\code{\link{OptPath}}].
#' @export
#' @examples
#' ps = makeParamSet(
#'   makeIntegerParam("u", trafo=function(x) 2*x),
#'   makeNumericVectorParam("v", len=2, trafo=function(x) x/sum(x)),
#'   makeDiscreteParam("w", values=c("a", "b"))
#' )
#' op = makeOptPathDF(ps, y.names="y", minimize=TRUE)
#' addOptPathEl(op, x=list(3, c(2, 4), "a"), y=0, dob=1, eol=1)
#' addOptPathEl(op, x=list(4, c(5, 3), "b"), y=2, dob=5, eol=7)
#'
#' as.data.frame(op)
#' op = trafoOptPath(op)
#' as.data.frame(op)
trafoOptPath = function(opt.path) {
  assertClass(opt.path, "OptPath")
  if (opt.path$add.transformed.x)
    stop("Cannot further trafo opt.path, you already added transformed x values to it!")
  ps = opt.path$par.set
  # FIXME: this only works for the DF implementation!
  op2 = makeOptPathDF(opt.path$par.set, opt.path$y.names, opt.path$minimize, add.transformed.x = TRUE)
  lapply(1:getOptPathLength(opt.path), function(i) {
    z = getOptPathEl(opt.path, i)
    x = trafoValue(ps, z$x)
    addOptPathEl(op2, x, z$y, z$dob, z$eol)
  })
  return(op2)
}
