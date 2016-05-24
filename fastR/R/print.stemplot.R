#' @export
print.stemplot <-
function (x, ...) 
{
    for (i in seq(x)) cat(x[[i]], sep = "\n")
    invisible(x)
}
