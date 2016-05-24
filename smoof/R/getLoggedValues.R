#' Extract logged values of a function wrapped by a logging wrapper.
#'
#' @param fn [\code{wrapped_smoof_function}]\cr
#'   Wrapped smoof function.
#' @param compact [\code{logical(1)}]\cr
#'   Wrap all logged values in a single data frame? Default is \code{FALSE}.
#' @return [\code{list} || \code{data.frame}]
#'   If \code{compact} is \code{TRUE}, a single data frame. Otherwise the function
#'   returns a list containing the following values:
#'   \describe{
#'     \item{pars}{Data frame of parameter values, i.e., x-values or the empty
#'     data frame if x-values were not logged.}
#'     \item{obj.vals}{Numeric vector of objective vals in the single-objective
#'     case respectively a matrix of objective vals for multi-objective functions.}
#'   }
#' @seealso \code{\link{addLoggingWrapper}}
#' @export
getLoggedValues = function(fn, compact = FALSE) {
  UseMethod("getLoggedValues")
}

#' @export
getLoggedValues.smoof_logging_function = function(fn, compact = FALSE) {
  env = environment(fn)
  pars = env$pars
  obj.vals = env$obj.vals
  # wrap everything up in a single data frame
  if (compact) {
    # if only the x-values are stored just return the data frame
    if (is.null(obj.vals)) {
      return(pars)
    }
    df = as.data.frame(t(obj.vals))
    names(df) = paste0("y", seq(ncol(df)))
    if (!is.null(pars)) {
      # append x-values if stored
      df = cbind(pars, df)
    }
    return(df)
  }
  if (!is.null(obj.vals) && nrow(obj.vals) == 1L) {
    obj.vals = as.numeric(obj.vals)
  }
  return(list(pars = pars, obj.vals = obj.vals))
}

#' @export
getLoggedValues.smoof_wrapped_function = function(fn, compact = FALSE) {
  return(getLoggedValues(getWrappedFunction(fn), compact))
}
