#' Add Comments to Character Strings.
#'
#' This operator allows you to add comments to character strings.
#' 
#' @param x A character string that is to be commented.
#' @param y A character string (the comment).
#' @return Returns a character string of the class \code{subcom} with a comment 
#' added as a \code{"comment"} attribute.
#' @keywords comment
#' @export
#' @note The operator, \code{\%:)\%}, is a simple smiley face emotion because 
#' commented code is happy code.
#' @seealso \code{\link[base]{comment}}
#' @rdname add_comments
#' @examples
#' a <- "The character string"
#' b <- "The comment"
#' 
#' (out <- a %:)% b)
#' attributes(out)
#' comment(out)
#' 
#' minimal <- construct("a", "b", "c" %:)% "A love note to your future self")
#' minimal
#' comments(minimal)
`%:)%` <- function(x, y) { 
    class(x) <- c("subcom", "character")
    attributes(x)[["comment"]] <- y
    x
}

#' @export
#' @rdname add_comments
`%comment%` <- `%:)%`
