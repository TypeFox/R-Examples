#' Improved list of objects
#' 
#' Improved list of objects.  Sorts by size by default.  This was shamelessly 
#' stolen from \url{http://stackoverflow.com/q/1358003/654296}.
#' 
#' @author Dirk Eddelbuettel, Tony Breyal
#' @keywords NA
#'   
#' @param pos numeric. Position in the stack.
#' @param pattern Regex to filter the objects by.
#' @param order.by character. Either 'Type', 'Size', 'PrettySize', 'Rows', or 
#'   'Columns'. This will dictate how the output is ordered.
#' @param decreasing logical. Should the output be displayed in decreasing order?
#' @param head logical. Use head on the output?
#' @param n numeric. Number of objects to display is head is TRUE.
#'   
#' @return A data.frame with type, size in bytes, human-readable size, rows, and
#'   columns of every object in the environment.
#' @import utils
#' 
#' @examples
#' \dontrun{
#' a <- rnorm(100000)
#' b <- matrix(1, 1000, 100)
#' lsa()
#' }
#'   
#' @export
lsa <- function (pos = 1, pattern, order.by = "Size",
                  decreasing=TRUE, head=TRUE, n=10) {
    napply <- function(names, fn) sapply(names, function(x)
        fn(get(x, pos = pos)))
    names <- ls(pos = pos, pattern = pattern)
    obj.class <- napply(names, function(x) as.character(class(x))[1])
    obj.mode <- napply(names, mode)
    obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
    obj.prettysize <- napply(names, function(x) {
        capture.output(print(object.size(x), units = "auto")) })
    obj.size <- napply(names, object.size)
    obj.dim <- t(napply(names, function(x)
        as.numeric(dim(x))[1:2]))
    vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
    obj.dim[vec, 1] <- napply(names, length)[vec]
    out <- data.frame(obj.type, obj.size, obj.prettysize, obj.dim)
    names(out) <- c("Type", "Size", "PrettySize", "Rows", "Columns")
    if (!missing(order.by))
        out <- out[order(out[[order.by]], decreasing=decreasing), ]
    if (head)
        out <- head(out, n)
    out
}