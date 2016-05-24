#' Get/Set Comments From a regexr Object
#' 
#' \code{comments} - Get the \code{comments} from a \code{regexr} object.
#' 
#' @param x A regexr object.
#' @param value The comment(s) to assign.
#' @param \ldots Ignored.
#' @rdname comments
#' @export
#' @return \code{comments} - Returns a list of comments.
#' @examples 
#' minimal <- construct("a", "b", "c" %:)% "Comment #3")
#' minimal
#' comments(minimal)
#' comments(minimal)[2] <- "A comment"
#' comments(minimal)
#' 
#' minimal <- construct("a", "b", "c")
#' out <- set_comments(minimal, paste("comment", 1:3))
#' comments(out) 
comments <- function (x, ...){
    UseMethod("comments")
}

#' Comments of a regexr Object
#' 
#' \code{comments<-} - Set the \code{comments} of a \code{regexr} object.
#' 
#' @rdname comments
#' @export
`comments<-` <- function(x, value){
    UseMethod("comments<-")
}


#' Set the Comments in a \code{regexr} Object
#' 
#' \code{set_comments} - This is a convenience function that sets the 
#' \code{\link[regexr]{comments}} on a \code{regexr} object and returns the 
#' object.
#' 
#' @param y The comments to assign.
#' @return \code{set_comments} - Returns a \code{regexr} object.
#' @export
#' @rdname comments
set_comments <- function (x, y) {
    comments(x) <- y
    x
}

