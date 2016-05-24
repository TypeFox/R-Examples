#' Load RData file and return objects in it.
#'
#' @param file [\code{character(1)}]\cr
#'   File to load.
#' @param parts [\code{character}]\cr
#'   Elements in file to load.
#'   Default is all.
#' @param simplify [\code{logical(1)}]\cr
#'   If \code{TRUE}, a list is only returned if \code{parts} and the file contain both more
#'   than 1 element, otherwise the element is directly returned.
#'   Default is \code{TRUE}.
#' @param envir [\code{environment(1)}]\cr
#'   Assign objects to this environment.
#'   Default is not to assign.
#' @param impute [\code{ANY}]\cr
#'   If \code{file} does not exists, return \code{impute} instead.
#'   Default is missing which will result in an exception if \code{file} is not found.
#' @return Either a single object or a list.
#' @export
#' @examples
#' fn = tempfile()
#' save2(file = fn, a = 1, b = 2, c = 3)
#' load2(fn, parts = "a")
#' load2(fn, parts = c("a", "c"))
load2 = function(file, parts, simplify = TRUE, envir, impute) {
  assertFlag(simplify)
  ee = new.env()
  if (!missing(impute) && !file.exists(file))
    return(impute)
  load(file, envir = ee)
  ns = ls(ee, all.names = TRUE)
  if (!missing(parts)) {
    assertCharacter(parts, any.missing = FALSE)
    d = setdiff(parts, ns)
    if (length(d) > 0L)
      stopf("File %s does not contain: %s", file, collapse(d))
  } else {
    parts = ns
  }
  if (!missing(envir)) {
    lapply(ns, function(x) assign(x, ee[[x]], envir = envir))
  }
  if (simplify) {
    if (length(ns) == 1L)
      return(ee[[ns]])
    if (length(parts) == 1L)
      return(ee[[parts]])
  }
  mget(parts, envir = ee)
}
