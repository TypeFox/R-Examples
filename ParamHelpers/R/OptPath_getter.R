#' Get the length of the optimization path.
#'
#' Dependent parameters whose requirements are not satisfied are represented by a scalar
#' NA in the output.
#'
#' @template arg_op
#' @return [\code{integer(1)}]
#' @export
#' @family optpath
getOptPathLength = function(op) {
  UseMethod("getOptPathLength")
}

#' Get an element from the optimization path.
#'
#' Dependent parameters whose requirements are not satisfied are represented by a scalar NA
#' in the elements of \code{x} of the return value.
#'
#' @template arg_op
#' @param index [\code{integer(1)}]\cr
#'   Index of element.
#' @return List with elements \code{x} [named \code{list}], \code{y} [named \code{numeric}],
#'   \code{dob} [\code{integer(1)}], \code{eol} [\code{integer(1)}].
#'   The elements \code{error.message} [\code{character(1)}],
#'   \code{exec.time} [\code{numeric(1)}] and \code{extra} [named \code{list}] are
#'   there if the respective options in \code{\link{OptPath}} are enabled.
#' @rdname getOptPathEl
#' @export
#' @family optpath
getOptPathEl = function(op, index) {
  UseMethod("getOptPathEl")
}

#' Get data.frame of input points (X-space) referring to the param set from the optimization path.
#'
#' @template arg_op
#' @template arg_opgetter_dob
#' @template arg_opgetter_eol
#' @return [\code{data.frame}].
#' @export
#' @family optpath
getOptPathX = function(op, dob, eol) {
  UseMethod("getOptPathX")
}

#' Get y-vector or y-matrix from the optimization path.
#'
#' @template arg_op
#' @param names [\code{character}]\cr
#'   Names of performance measure.
#'   Default is all performance measures in path.
#' @template arg_opgetter_dob
#' @template arg_opgetter_eol
#' @param drop [\code{logical(1)}]\cr
#'   Return vector instead of matrix when only one y-column was selected?
#'   Default is \code{TRUE}.
#' @return [\code{numeric} | \code{matrix}]. The columns of the matrix are always named.
#' @export
#' @family optpath
getOptPathY = function(op, names, dob, eol, drop = TRUE) {
  UseMethod("getOptPathY")
}

#' Get date-of-birth vector from the optimization path.
#'
#' @template arg_op
#' @template arg_opgetter_dob
#' @template arg_opgetter_eol
#' @return [\code{integer}].
#' @export
#' @family optpath
getOptPathDOB = function(op, dob, eol) {
  UseMethod("getOptPathDOB")
}

#' Get end-of-life vector from the optimization path.
#'
#' @template arg_op
#' @template arg_opgetter_dob
#' @template arg_opgetter_eol
#' @return [\code{integer}].
#' @export
#' @family optpath
getOptPathEOL = function(op, dob, eol) {
  UseMethod("getOptPathEOL")
}

#' Get error-message vector from the optimization path.
#'
#' @template arg_op
#' @template arg_opgetter_dob
#' @template arg_opgetter_eol
#' @return [\code{character}].
#' @export
#' @family optpath
getOptPathErrorMessages = function(op, dob, eol) {
  UseMethod("getOptPathErrorMessages")
}

#' Get exec-time vector from the optimization path.
#'
#' @template arg_op
#' @template arg_opgetter_dob
#' @template arg_opgetter_eol
#' @return [\code{numeric}].
#' @export
#' @family optpath
getOptPathExecTimes = function(op, dob, eol) {
  UseMethod("getOptPathExecTimes")
}

#' Get column from the optimization path.
#'
#' @template arg_op
#' @param name [\code{character(1)}]\cr
#'   Name of the column.
#' @template arg_opgetter_dob
#' @template arg_opgetter_eol
#' @return Single column as a vector.
#' @export
#' @family optpath
getOptPathCol = function(op, name, dob, eol) {
  UseMethod("getOptPathCol")
}

#' Get columns from the optimization path.
#'
#' @template arg_op
#' @param names [\code{character}]\cr
#'   Names of the columns.
#' @template arg_opgetter_dob
#' @template arg_opgetter_eol
#' @inheritParams as.data.frame.OptPathDF
#' @return [\code{data.frame}].
#' @export
#' @family optpath
getOptPathCols = function(op, names, dob, eol, row.names = NULL) {
  UseMethod("getOptPathCols")
}

#' Get index of the best element from optimization path.
#'
#' @template arg_op
#' @param y.name [\code{character(1)}]\cr
#'   Name of target value to decide which element is best.
#'   Default is \code{y.names[1]}.
#' @template arg_opgetter_dob
#' @template arg_opgetter_eol
#' @param ties [\code{character(1)}]\cr
#'   How should ties be broken when more than one optimal element is found?
#'   \dQuote{all}: return all indices,
#'   \dQuote{first}: return first optimal element in path,
#'   \dQuote{last}: return last optimal element in path,
#'   \dQuote{random}: return random optimal element in path.
#'   Default is \dQuote{last}.
#' @return [\code{integer}]
#'   Index or indices into path. See \code{ties}.
#' @export
#' @family optpath
#' @examples
#' ps = makeParamSet(makeNumericParam("x"))
#' op = makeOptPathDF(par.set = ps, y.names = "y", minimize = TRUE)
#' addOptPathEl(op, x = list(x = 1), y = 5)
#' addOptPathEl(op, x = list(x = 2), y = 3)
#' addOptPathEl(op, x = list(x = 3), y = 9)
#' addOptPathEl(op, x = list(x = 4), y = 3)
#' as.data.frame(op)
#' getOptPathBestIndex(op)
#' getOptPathBestIndex(op, ties = "first")
getOptPathBestIndex = function(op, y.name = op$y.names[1], dob = op$env$dob, eol = op$env$eol, ties = "last") {
  assertClass(op, "OptPath")
  assertChoice(y.name, choices = op$y.names)
  dob = asInteger(dob, any.missing = TRUE)
  eol = asInteger(eol, any.missing = TRUE)
  assertChoice(ties, c("all", "first", "last", "random"))
  life.inds = which(op$env$dob %in% dob & op$env$eol %in% eol)
  if (length(life.inds) == 0)
    stop("No element found which matches dob and eol restrictions!")
  y = getOptPathY(op, y.name)[life.inds]
  if (all(is.na(y))) {
    best.inds = life.inds
  } else {
    if (op$minimize[y.name])
      best.inds = which(min(y, na.rm = TRUE) == y)
    else
      best.inds = which(max(y, na.rm = TRUE) == y)
    best.inds = life.inds[best.inds]
  }
  if (length(best.inds) > 1) {
    if (ties == "all")
      return(best.inds)
    else if (ties == "first")
      return(best.inds[1])
    else if (ties == "last")
      return(best.inds[length(best.inds)])
    else if (ties == "random")
      return(best.inds[sample(length(best.inds), 1)])
  } else {
    return(best.inds)
  }
}

#' Get indices of pareto front of optimization path.
#'
#' @template arg_op
#' @param y.names [\code{character}]\cr
#'   Names of performance measures to construct pareto front for.
#'   Default is all performance measures.
#' @template arg_opgetter_dob
#' @template arg_opgetter_eol
#' @param index [\code{logical(1)}]\cr
#'   Return indices into path of front or y-matrix of nondominated points?
#'   Default is \code{FALSE}.
#' @return [\code{matrix} | \code{integer}]. Either matrix (with named columns) of points of front
#'   in objective space or indices into path for front.
#' @export
#' @family optpath
#' @examples
#' ps = makeParamSet(makeNumericParam("x"))
#' op = makeOptPathDF(par.set = ps, y.names = c("y1", "y2"), minimize = c(TRUE, TRUE))
#' addOptPathEl(op, x = list(x = 1), y = c(5, 3))
#' addOptPathEl(op, x = list(x = 2), y = c(2, 4))
#' addOptPathEl(op, x = list(x = 3), y = c(9, 4))
#' addOptPathEl(op, x = list(x = 4), y = c(4, 9))
#' as.data.frame(op)
#' getOptPathParetoFront(op)
#' getOptPathParetoFront(op, index = TRUE)
getOptPathParetoFront = function(op, y.names = op$y.names, dob = op$env$dob, eol = op$env$eol, index = FALSE) {
  assertClass(op, "OptPath")
  assertCharacter(y.names, min.len = 2)
  assertSubset(y.names, op$y.names, empty.ok = FALSE)
  dob = asInteger(dob, any.missing = TRUE)
  eol = asInteger(eol, any.missing = TRUE)
  assertFlag(index, na.ok = TRUE)
  requirePackages("emoa", default.method = "load")
  life.inds = which(op$env$dob %in% dob & op$env$eol %in% eol)
  if (length(life.inds) == 0)
    stop("No element found which matches dob and eol restrictions!")
  y = getOptPathY(op, y.names, drop = FALSE)[life.inds, , drop = FALSE]
  # multiply columns with -1 if maximize
  k = ifelse(op$minimize, 1, -1)
  y2 = t(y) * k
  # is_dominated has kind of buggy behavoiur if y2 is a row
  # (it hinks, we have a 1-dimensional optimization prob und returns the min index)
  # so we have to treat this case manually
  if (nrow(y2) == 1)
    nondom = 1
  else
    nondom = which(!emoa::is_dominated(y2))
  if (index)
    return(life.inds[nondom])
  else
    return(y[nondom, , drop = FALSE])
}

