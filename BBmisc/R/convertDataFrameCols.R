#' Converts columns in a data frame to characters, factors or numerics.
#'
#' @param df [\code{data.frame}]\cr
#'   Data frame.
#' @param chars.as.factor [\code{logical(1)}]\cr
#'   Should characters be converted to factors?
#'   Default is \code{FALSE}.
#' @param factors.as.char [\code{logical(1)}]\cr
#'   Should characters be converted to factors?
#'   Default is \code{FALSE}.
#' @param ints.as.num [\code{logical(1)}]\cr
#'   Should integers be converted to numerics?
#'   Default is \code{FALSE}.
#' @param logicals.as.factor [\code{logical(1)}]\cr
#'   Should logicals be converted to factors?
#'   Default is \code{FALSE}.
#' @export
#' @return [\code{data.frame}].
convertDataFrameCols = function(df, chars.as.factor = FALSE, factors.as.char = FALSE, ints.as.num = FALSE, logicals.as.factor = FALSE) {
  assertDataFrame(df)
  assertFlag(chars.as.factor)
  assertFlag(factors.as.char)
  assertFlag(ints.as.num)
  assertFlag(logicals.as.factor)
  df = x = as.list(df)

  if (chars.as.factor) {
    i = vlapply(df, is.character)
    if (any(i))
      x[i] = lapply(x[i], factor)
  }

  if (factors.as.char) {
    i = vlapply(df, is.factor)
    if (any(i))
      x[i] = lapply(x[i], as.character)
  }

  if (ints.as.num) {
    i = vlapply(df, is.integer)
    if (any(i))
      x[i] = lapply(x[i], as.double)
  }

  if (logicals.as.factor) {
    i = vlapply(df, is.logical)
    if (any(i))
      x[i] = lapply(x[i], factor, levels = c("TRUE", "FALSE"))
  }

  as.data.frame(x, stringsAsFactors = FALSE)
}
