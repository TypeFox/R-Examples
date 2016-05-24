#' Create a source abstraction interface layer (SAIL) object.
#'
#' This function returns an object of class \code{sail} which behaves
#' like \code{\link{fail}}, but is indented for loading and saving
#' R source code files.
#'
#' @param path [\code{character(1)}]\cr
#'   Path to work in, will be created if it does not exists.
#' @param extension [\code{character(1)}]\cr
#'   File extension to work with.
#'   Default is \dQuote{R}.
#' @param all.files [\code{logical(1)}]\cr
#'   Also include hidden files, i.e. files whose name start with a dot (\dQuote{.}).
#'   Default is \code{FALSE}.
#' @param use.cache [\code{logical(1)}]\cr
#'   Use a memory cache per global default.
#'   Global option which can locally be overwritten in most functions.
#'   Default is \code{FALSE}
#' @param simplify [\code{character(1)}]\cr
#'   If only one object is defined in a sourced R file,
#'   should the return value be simplified? If set to \code{TRUE},
#'   instead of a list containing one element the object itself will be returned.
#' @param suppressMessages [\code{logical(1)}]\cr
#'   Wrap the \code{\link[base]{sys.source}} command into \code{\link[base]{suppressMessages}}
#'   and \code{link[base]{suppressPackageStartupMessages}}?
#'   Default is \code{FALSE}, i.e. you will see regular output of sourced scripts.
#' @return Object of class \code{sail}. See the documentation of \code{\link{fail}}
#'   for details.
#' @export
sail = function(path = getwd(), extension = "R", all.files = FALSE, use.cache = FALSE, simplify = TRUE,
  suppressMessages = FALSE) {
  .self = list(
    path = checkPath(path),
    extension = checkExtension(extension),
    all.files = asFlag(all.files),
    use.cache = asFlag(use.cache),
    simplify = asFlag(simplify, na.ok = TRUE),
    cache = Cache(),
    loadFun = loadR,
    saveFun = saveR,
    suppressMessages = asFlag(suppressMessages)
  )
  checkCollision(Ls(.self))
  setClasses(makeObject(.self), "sail")
}
