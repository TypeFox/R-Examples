find.param <- function (object) {
    if (!inherits(object, 'estimatetables'))
        stop("input unsuitable")
    for (i in 1: length(object$output)) {
        tmp <- row.names(object$output[[1]][[1]])
        if (is.character(tmp))
            break
    }
    tmp
}

find.stats <- function (object) {
    if (!inherits(object, 'selectedstatistics'))
        stop("input unsuitable")
    for (i in 1: length(object$output)) {
        tmp <- dimnames(object$output[[1]])[[2]]
        if (is.character(tmp))
            break
    }
    tmp
}

