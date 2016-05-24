#' Get Sub-expressions from \code{regexr} Object
#' 
#' Get sub-expressions from \code{regexr} object.
#' 
#' @param x A \code{regexr} object.
#' @param \ldots Ignored.
#' @export
#' @return Returns a list of regular expression chunks.
#' @examples 
#' minimal <- construct("a", "b", "c")
#' minimal
#' unglue(minimal)
unglue <- function (x, ...){
    UseMethod("unglue")
}