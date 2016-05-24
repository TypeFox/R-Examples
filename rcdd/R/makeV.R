makeV <- function(points, rays, lines, x = NULL) {

    if (missing(points) && missing(rays) && missing(lines))
       stop("at least one of 'points' and 'rays' and 'lines' must be specified")

    rational.output <- (! missing(points)) && is.character(points) ||
        (! missing(rays)) && is.character(rays) ||
        (! missing(lines)) && is.character(lines) ||
        (! is.null(x)) && is.character(x)

    d <- 0
    n <- 0
    if (! missing(points)) {
        stopifnot(is.numeric(points) || is.character(points))
        if (! is.matrix(points))
            points <- rbind(points)
        d <- ncol(points)
        n <- nrow(points)
        if (is.numeric(points)) {
            stopifnot(is.finite(points))
            if (rational.output)
                points <- d2q(points)
        } else {
            fubar <- try(q2q(points))
            if(inherits(fubar, "try-error"))
                stop("'points' character but not GMP rational")
        }
    }
    if (! missing(rays)) {
        stopifnot(is.numeric(rays) || is.character(rays))
        if (! is.matrix(rays))
            rays <- rbind(rays)
        if (d > 0) {
            if (d != ncol(rays))
                stop("column dimensions of arguments differ")
        } else {
            d <- ncol(rays)
        }
        n <- n + nrow(rays)
        if (is.numeric(rays)) {
            stopifnot(is.finite(rays))
            if (rational.output)
                rays <- d2q(rays)
        } else {
            fubar <- try(q2q(rays))
            if(inherits(fubar, "try-error"))
                stop("'rays' character but not GMP rational")
        }
    }
    if (! missing(lines)) {
        stopifnot(is.numeric(lines) || is.character(lines))
        if (! is.matrix(lines))
            lines <- rbind(lines)
        if (d > 0) {
            if (d != ncol(lines))
                stop("column dimensions of arguments differ")
        } else {
            d <- ncol(lines)
        }
        n <- n + nrow(lines)
        if (is.numeric(lines)) {
            stopifnot(is.finite(lines))
            if (rational.output)
                lines <- d2q(lines)
        } else {
            fubar <- try(q2q(lines))
            if(inherits(fubar, "try-error"))
                stop("'lines' character but not GMP rational")
        }
    }
    if (! is.null(x)) {
        stopifnot(is.numeric(x) || is.character(x))
        stopifnot(is.matrix(x))
        if (d > 0 && d + 2 != ncol(x))
            stop("column dimension of 'x' not 2 + column dimensions of other arguments")
        if (is.numeric(x)) {
            stopifnot(is.finite(x))
            if (rational.output)
                x <- d2q(x)
        } else {
            fubar <- try(q2q(x))
            if(inherits(fubar, "try-error"))
                stop("'x' character but not GMP rational")
        }
        n <- n + nrow(x)
    }
    if (n == 0) stop("all arguments have row dimension zero")

    fred <- attr(x, "representation")
    if ((! is.null(fred)) && (fred != "V"))
        stop("\"representation\" attribute of argument 'x' not \"V\"")

    foo <- NULL
    if (! missing(points))
        foo <- cbind(0, 1, points, deparse.level = 0)

    bar <- NULL
    if (! missing(rays))
        bar <- cbind(0, 0, rays, deparse.level = 0)

    baz <- NULL
    if (! missing(lines))
        baz <- cbind(1, 0, lines, deparse.level = 0)

    qux <- rbind(x, foo, bar, baz, deparse.level = 0)
    attr(qux, "representation") <- "V"
    validcdd(qux)
    return(qux)
}

addVpoints <- function(points, x) makeV(points = points, x = x)

addVrays <- function(rays, x) makeV(rays = rays, x = x)

addVlines <- function(lines, x) makeV(lines = lines, x = x)

