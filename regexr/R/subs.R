#' Get/Set Regex Sub-expressions From a regexr Object
#' 
#' \code{subs} - Get the sub-expressions from a \code{regexr} object.
#' 
#' @param x A regexr object.
#' @param value The comment(s) to assign.
#' @param \ldots Ignored.
#' @rdname subs
#' @export
#' @return \code{subs} - Returns a list of sub-expressions.
#' @examples 
#' minimal <- construct("a", "b", "c")
#' minimal
#' subs(minimal)
#' subs(minimal)[2] <- "\\s+[A-Z]|[0-9]"
#' subs(minimal)
#' 
#' minimal <- construct("a", "b", "c")
#' out <- set_subs(minimal, c("(", "\\s??", ")"))
#' subs(out)
subs <- function (x, ...){
    UseMethod("subs")
}

#' Set Regex Sub-expressions From a regexr Object
#' 
#' \code{subs<-} - Set the sub-expressions(s) of a \code{regexr} object.
#' 
#' @rdname subs
#' @export
`subs<-` <- function(x, value){
    UseMethod("subs<-")
}


#' Set the Sub-expressions in a \code{regexr} Object
#' 
#' \code{set_subs} - This is a convenience function that sets the 
#' \code{\link[regexr]{subs}} on a \code{regexr} object and returns the object. 
#' 
#' @param y The sub-expressions to assign.
#' @return \code{set_subs} - Returns a \code{regexr} object.
#' @export
#' @rdname subs
set_subs <- function (x, y) {
    subs(x) <- y
    x
}

