compose <- function(x, 
                    y,
                    e=NULL,
                    alg=c('goedel', 'goguen', 'lukasiewicz'),
                    type=c('basic', 'sub', 'super', 'square'),
                    quantifier=NULL) {

    stopifnot(is.matrix(x))
    stopifnot(is.matrix(y))
    stopifnot(is.numeric(x))
    stopifnot(is.numeric(y))
    stopifnot(nrow(x) > 0)
    stopifnot(ncol(y) > 0)
    stopifnot(ncol(x) == nrow(y))

    if (!is.null(e)) {
        stopifnot(is.matrix(e))
        stopifnot(is.numeric(e))
        stopifnot(nrow(y) == nrow(e))
        stopifnot(ncol(y) == ncol(e))
    }

    if (is.character(alg)) {
        alg <- match.arg(alg)
        alg <- algebra(alg)
    }
    stopifnot(is.algebra(alg))

    if (is.character(type)) {
        type <- match.arg(type)
        if (type == 'basic') {
            type <- alg$pt
            merge <- goedel.tconorm
        } else if (type == 'sub') {
            type <- alg$r
            merge <- goedel.tnorm
        } else if (type == 'super') {
            type <- function(x, y) { alg$r(y, x) }
            merge <- goedel.tnorm
        } else if (type == 'square') {
            type <- alg$b
            merge <- goedel.tnorm
        } else {
            stop('Unrecognized composition type')
        }
    }
    stopifnot(is.function(type))

    if (is.function(quantifier)) {
        merge <- function(x) {
            res <- sort(x, decreasing=TRUE)
            relcard <- seq_along(res) / length(res)
            goedel.tconorm(alg$pt(res, quantifier(relcard)))
        }
    } else if (!is.null(quantifier)) {
        stop("'quantifier' must be a function or NULL")
    }

    f <- function(x, y) {
        merge(type(x, y))
    }

    res <- mult(x, y, f)
    if (!is.null(e)) {
        warning('Computations with an excluding relation "e" is experimental!')
        e <- mult(x, e, f)
        res <- alg$pt(res, alg$n(e))
    }
    return(res)
}
