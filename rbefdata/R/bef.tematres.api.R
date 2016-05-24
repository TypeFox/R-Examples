#' Basic tematres API via rtematres
#'
#' This function queries a tematres thesaurus via the rtematres api package.
#'
#' @param ... Any task or term to be passed to the rtematres package
#' @import rtematres
#' @export bef.tematres.api

bef.tematres.api <- function(...) {
      rtematres.api.do(...)
}
