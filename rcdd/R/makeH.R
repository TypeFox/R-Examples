makeH <- function(a1, b1, a2, b2, x = NULL) {

    if (missing(a1) != missing(b1))
        stop("must have either both or neither of 'a1' and 'b1' missing")
    if (missing(a2) != missing(b2))
        stop("must have either both or neither of 'a2' and 'b2' missing")
    if (missing(a1) && missing(a2))
        stop("'a1' and 'b1' and 'a2' and 'b2' all missing")
    if ((! missing(a1)) && (! is.matrix(a1)))
        a1 <- matrix(a1, nrow = 1)
    if ((! missing(a2)) && (! is.matrix(a2)))
        a2 <- matrix(a2, nrow = 1)
    if ((! missing(a1)) && (! missing(a2)))
        if (ncol(a1) != ncol(a2))
            stop("ncol(a1) != ncol(a2)")
    if (! missing(a1)) {
        stopifnot(nrow(a1) == length(b1))
        stopifnot(is.numeric(a1) || is.character(a1))
        stopifnot(is.numeric(b1) || is.character(b1))
        if (is.numeric(a1)) {
            stopifnot(is.finite(a1))
        } else {
            fubar <- try(q2q(a1))
            if(inherits(fubar, "try-error"))
                stop("'a1' character but not GMP rational")
        }
        if (is.numeric(b1)) {
            stopifnot(is.finite(b1))
        } else {
            fubar <- try(q2q(b1))
            if(inherits(fubar, "try-error"))
                stop("'b1' character but not GMP rational")
        }
    }
    if (! missing(a2)) {
        stopifnot(nrow(a2) == length(b2))
        stopifnot(is.numeric(a2) || is.character(a2))
        stopifnot(is.numeric(b2) || is.character(b2))
        if (is.numeric(a2)) {
            stopifnot(is.finite(a2))
        } else {
            fubar <- try(q2q(a2))
            if(inherits(fubar, "try-error"))
                stop("'a2' character but not GMP rational")
        }
        if (is.numeric(b2)) {
            stopifnot(is.finite(b2))
        } else {
            fubar <- try(q2q(b2))
            if(inherits(fubar, "try-error"))
                stop("'b2' character but not GMP rational")
        }
    }

    fred <- attr(x, "representation")
    if ((! is.null(fred)) && (fred != "H"))
        stop("\"representation\" attribute of argument 'x' not \"H\"")
    if ((! is.null(x)) && (! is.matrix(x)))
        stop("argument 'x' must be NULL or matrix")
    if (! is.null(x)) {
        stopifnot(is.numeric(x) || is.character(x))
        if (is.numeric(x)) {
            stopifnot(is.finite(x))
        } else {
            fubar <- try(q2q(x))
            if(inherits(fubar, "try-error"))
                stop("'x' character but not GMP rational")
        }
    }

    rational.output <-
        ((! missing(a1)) && (is.character(a1) || is.character(b1))) ||
        ((! missing(a2)) && (is.character(a2) || is.character(b2))) ||
        ((! is.null(x)) && is.character(x))

    if (! rational.output) {

    foo <- NULL
    if (! missing(a1))
        foo <- cbind(0, as.vector(b1), - a1)

    bar <- NULL
    if (! missing(a2))
        bar <- cbind(1, as.vector(b2), - a2)

    } else {

        foo <- NULL
        if (! missing(a1)) {
            if (is.numeric(a1)) a1 <- d2q(a1)
            if (is.numeric(b1)) b1 <- d2q(b1)
            foo <- cbind(0, as.vector(b1), qneg(a1))
        }
        bar <- NULL
        if (! missing(a2)) {
            if (is.numeric(a2)) a2 <- d2q(a2)
            if (is.numeric(b2)) b2 <- d2q(b2)
            bar <- cbind(1, as.vector(b2), qneg(a2))
        }
        if (! is.null(x)) {
            if (is.numeric(x)) x <- d2q(x)
        }
    }

    mcol <- NULL
    if (! is.null(x))
        mcol <- ncol(x)
    if (! is.null(foo)) {
        if (is.null(mcol)) {
            mcol <- ncol(foo)
        } else {
            if (mcol != ncol(foo))
                stop("ncol(a1) + 2 != ncol(x)")
        }
    }
    if (! is.null(bar)) {
        if (is.null(mcol)) {
            mcol <- ncol(bar)
        } else {
            if (mcol != ncol(bar))
                stop("ncol(a2) + 2 != ncol(x) || ncol(a2) != ncol(a1)")
        }
    }
    baz <- rbind(x, bar, foo)
    dimnames(baz) <- NULL
    attr(baz, "representation") <- "H"
    validcdd(baz)
    return(baz)
}

addHeq <- function(a, b, x) makeH(a2 = a, b2 = b, x = x)

addHin <- function(a, b, x) makeH(a1 = a, b1 = b, x = x)

