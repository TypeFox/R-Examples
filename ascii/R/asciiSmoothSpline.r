##' @export
##' @method ascii smooth.spline
ascii.smooth.spline <- function(x, caption = NULL, caption.level = NULL, list.type = "bullet", ...) {
    x <- as.list(capture.output(x)[-1:-3])
    obj <- asciiList$new(x = x, caption = caption, caption.level = caption.level, list.type = list.type)
    return(obj)
}

