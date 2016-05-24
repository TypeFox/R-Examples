namer <- function(x, ...){
    if (is.null(names(x))) names(x) <- rep("", length(x))
    x
}
 

get_comment <- function(x, ...) {
     attributes(x)[["comment"]]
}


is.regex <- function (pattern) {
    out <- suppressWarnings(try(gsub(pattern, "", "hello", perl = TRUE), 
        silent = TRUE))
    ifelse(inherits(out, "try-error"), FALSE, TRUE)
}










