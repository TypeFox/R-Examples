#' Set the Names in a \code{regexr} Object
#' 
#' This is a convenience function that sets the names on a \code{regexr} object 
#' and returns the object. This function works the same as 
#' \code{\link[stats]{setNames}} but provides a naming which is consistent with
#' \code{set_regex} and \code{set_comments}.
#' 
#' @param x The \code{regexr} object.
#' @param y The names to assign.
#' @return Returns a \code{regexr} object.
#' @export
#' @seealso \code{\link[stats]{setNames}}
#' @examples
#' minimal <- construct("a", "b", "c")
#' out <- set_names(minimal, 1:3)
#' names(out)
set_names <- function(x, y){
    names(x) <- y
    x
}
