#' Parses \code{...} arguments to a named list.
#'
#' The deparsed name will be used for arguments with missing names.
#' Missing names will be set to \code{NA}.
#'
#' @param ...
#'   Arbitrary number of objects.
#' @return [\code{list}]: Named list with objects.
#' @export
#' @examples
#' z = 3
#' argsAsNamedList(x = 1, y = 2, z)
argsAsNamedList = function(...) {
  args = list(...)
  ns = names2(args)
  ns.missing = is.na(ns)
  if (any(ns.missing)) {
    ns.sub = as.character(substitute(deparse(...)))[-1L]
    ns[ns.missing] = ns.sub[ns.missing]
  }
  setNames(args, replace(ns, ns %in% c("NA", "NULL", ""), NA_character_))
}
