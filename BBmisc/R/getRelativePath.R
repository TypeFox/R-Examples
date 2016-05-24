#' Construct a path relative to another
#'
#' Constructs a relative path from path \code{from} to path \code{to}.
#' If this is not possible (i.e. different drive letters on windows systems),
#' \code{NA} is returned.
#'
#' @param to [\code{character(1)}]\cr
#'  Where the relative path should point to.
#' @param from [\code{character(1)}]\cr
#'  From which part to start.
#'  Default is \code{\link[base]{getwd}}.
#' @param ignore.case [\code{logical(1)}]\cr
#'  Should path comparisons be made case insensitve?
#'  Default is \code{TRUE} on Windows systems and \code{FALSE} on other systems.
#' @return [character(1)]: A relative path.
#' @export
getRelativePath = function(to, from = getwd(), ignore.case = isWindows()) {
  numberCommonParts = function(p1, p2) {
    for (i in seq_len(min(length(p1), length(p2)))) {
      if (p1[i] != p2[i])
        return(i - 1L)
    }
    return(if (is.null(i)) 0L else i)
  }

  from = splitPath(from)
  to = splitPath(to)
  assertFlag(ignore.case)
  if (length(from$drive) != length(to$drive))
    return(NA_character_)
  if (length(from$drive) > 0L && length(to$drive) > 0L && from$drive != to$drive)
    return(NA_character_)

  if (ignore.case)
    i = numberCommonParts(tolower(from$path), tolower(to$path))
  else
    i = numberCommonParts(from$path, to$path)

  res = c(rep.int("..", length(from$path) - i), tail(to$path, ifelse(i == 0L, Inf, -i)))
  if (length(res) == 0L)
    res = "."
  collapse(res, .Platform$file.sep)
}
