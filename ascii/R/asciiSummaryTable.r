##' @export
##' @method ascii summary.table
ascii.summary.table <- function(x, caption = NULL, caption.level = NULL, list.type = "bullet", ...) {
    x <- as.list(capture.output(x))
    obj <- asciiList$new(x = x, caption = caption, caption.level = caption.level, list.type = list.type)
    return(obj)
}
