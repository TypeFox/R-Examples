`coca` <- function(y, ...) {
    if(is.null(class(y))) class(y) <- data.class(y)
    UseMethod("coca", y)
}
