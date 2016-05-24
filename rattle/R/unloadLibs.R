#' Unload pacakges
#'
#' Detach the list of pacakges, only detaching those that are on the
#' search path.
#'
#' @param l Vector of package names.
#' @return nothing.
#' @rdname unloadLibs
unloadLibs <- function(l)
{
  for (p in l)
  {
    pn <- sprintf("package:%s", p)
    if (pn %in% search()) detach(pn, character.only=TRUE)
  }
  invisible()
}
  
