#' @title Convert a data.frame row to list of parameter-value-lists.
#'
#' @description
#' Please note that (naturally) the columns of \code{df} have to be of the correct
#' type w.r.t. the corresponding parameter. The only exception are integer parameters
#' where the corresponding columns in \code{df} are allowed to be numerics.
#' And also see the argument \code{enforce.col.types} as a way around this restriction.
#'
#' \tabular{ll}{
#'  numeric(vector)   \tab  \code{numeric}  \cr
#'  integer(vector)   \tab  \code{integer}  \cr
#'  discrete(vector)  \tab  \code{factor} (names of values = levels) \cr
#'  logical(vector)   \tab  \code{logical}
#' }
#'
#' Dependent parameters whose requirements are not satisfied are represented by a scalar
#' NA in the output.
#'
#' @param df [\code{data.frame}]\cr
#'   Data.frame, potentially from \code{\link{OptPathDF}}.
#'   Columns are assumed to be in the same order as par.set.
#' @template arg_parset
#' @param i [\code{integer(1)}]\cr
#'   Row index.
#' @param enforce.col.types [\code{logical(1)}]\cr
#'   Should all \code{df} columns be initially converted to the type
#'   returned by \code{getParamTypes(df, df.cols = TRUE)}.
#'   This can help to work with \dQuote{non-standard} data.frames where the types are
#'   slightly \dQuote{off}. But note that there is no guarantee that this will
#'   work if the types are really wrong and there is no naturally correct way
#'   to convert them.
#'   Default is \code{FALSE}.
#' @param ... [any]\cr
#'   Arguments passed to \code{\link[BBmisc]{convertDataFrameCols}}
#' @return [\code{list}]. Named by parameter ids.
#' @export
#' @useDynLib ParamHelpers c_dfRowsToList
#' @rdname dfRowsToList
dfRowsToList = function(df, par.set, enforce.col.types = FALSE, ...) {
  assertClass(df, "data.frame")
  assertClass(par.set, "ParamSet")
  assertFlag(enforce.col.types)

  if (enforce.col.types) {
    types = getParamTypes(par.set, df.cols = TRUE)
    for (i in 1:length(types)) {
      tt = types[i]
      if (tt == "factor")
        df[,i] = as.factor(df[, i])
      else
        df[,i] = as(df[,i], tt)
    }
  }

  lens = getParamLengths(par.set)
  cnames = extractSubList(par.set$pars, "cnames", simplify = FALSE)
  int.type = convertTypesToCInts(getParamTypes(par.set, df.cols = TRUE))

  # factors to chars, so we can evaluate requires
  df = convertDataFrameCols(df, factors.as.char = TRUE, ...)
  # ints might just be encoded as nums in df, convert before going to C
  ints.as.double = mapply(function(type, col) type == 2L && is.double(col), type = int.type, col = df)
  df[ints.as.double] = lapply(df[ints.as.double], as.integer)

  .Call("c_dfRowsToList", df, par.set$pars, int.type, names(par.set$pars), lens, cnames, PACKAGE = "ParamHelpers")
}

#' @export
#' @rdname dfRowsToList
dfRowToList = function(df, par.set, i, enforce.col.types = FALSE, ...) {
  dfRowsToList(df = df, par.set = par.set, enforce.col.types = enforce.col.types, ...)[[i]]
}
