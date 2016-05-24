#' Prints a subcom object
#' 
#' Prints a subcom object.
#' 
#' @param x The subcom object.
#' @param \ldots Ignored.
#' @export
#' @method print subcom
print.subcom <- function(x, ...){
    class(x) <- "character"
    print(x)
}