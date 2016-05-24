
fam.bernoulli <- function() {
    result <- list(name = "bernoulli")
    class(result) <- "astfam"
    return(result)
}

fam.poisson <- function() {
    result <- list(name = "poisson")
    class(result) <- "astfam"
    return(result)
}

fam.zero.truncated.poisson <- function() {
    result <- list(name = "zero.truncated.poisson")
    class(result) <- "astfam"
    return(result)
}

fam.normal.location.scale <- function() {
    result <- list(name = "normal.location.scale")
    class(result) <- "astfam"
    return(result)
}

fam.multinomial <- function(dimension) {
    result <- list(name = "multinomial", dimension = dimension)
    class(result) <- "astfam"
    return(result)
}

as.character.astfam <- function(x, ...) {
    nam <- x$name
    if (is.null(nam))
        stop("astfam object with no name")
    x$name <- NULL

    if (length(x) == 0)
        return(nam)

    fred <- ""
    for (i in 1:length(x)) {
        if (fred != "")
            fred <- paste(fred, ", ", sep = "")
        fred <- paste(fred, names(x)[i], " = ", x[[i]], sep = "")
    }
    nam <- paste(nam, "(", fred, ")", sep = "")
    return(nam)
}

print.astfam <- function(x, ...) {
    foo <- as.character(x)
    cat(foo, "\n")
    return(invisible(foo))
}

fam.clear <- function() {
    .C("astfam_clear", PACKAGE = "aster2")
    return(invisible(NULL))
}

fam.set <- function(fam) {
    stopifnot(inherits(fam, "astfam"))
    name <- as.character(fam$name)
    fam$name <- NULL
    foo <- as.double(unlist(fam))
    if (! all(is.finite(foo)))
        stop("some hyperparameters not finite")
    if (length(foo) > 2)
        stop("more than 2 hyperparameters not (currently) allowed")
    hyper1 <- as.double(0)
    hyper2 <- as.double(0)
    if (length(foo) >= 1) hyper1 <- as.double(foo[[1]])
    if (length(foo) >= 2) hyper2 <- as.double(foo[[2]])
    .C("astfam_set", name, hyper1, hyper2, PACKAGE = "aster2")
    return(invisible(NULL))
}

fam.set.tolerance <- function(tolerance) {
    .C("astfam_set_tolerance", tolerance = as.double(tolerance),
        PACKAGE = "aster2")
    return(invisible(NULL))
}

fam.reset.tolerance <- function() {
    .C("astfam_reset_tolerance", PACKAGE = "aster2")
    return(invisible(NULL))
}

fam.dimension <- function(i) {
    stopifnot(is.atomic(i))
    stopifnot(is.numeric(i))
    stopifnot(i == as.integer(i))
    .C("astfam_dimension", fam = as.integer(i), dimen = integer(1),
        PACKAGE = "aster2")$dimen
}

cumulant <- function(theta, fam, deriv = 0, delta) {
    stopifnot(inherits(fam, "astfam"))
    stopifnot(is.atomic(theta))
    stopifnot(is.numeric(theta))
    stopifnot(all(is.finite(theta)))
    stopifnot(is.atomic(deriv))
    stopifnot(is.numeric(deriv))
    stopifnot(length(deriv) == 1)
    stopifnot(deriv == as.integer(deriv))
    stopifnot(deriv >= 0 && deriv <= 3)
    fam.clear()
    fam.set(fam)
    d <- fam.dimension(1)
    if (missing(delta)) delta <- rep(0, d)
    stopifnot(is.atomic(delta))
    stopifnot(is.numeric(delta))
    stopifnot(all(is.finite(delta)))
    if (length(theta) != d) stop("theta wrong dimension")
    if (length(delta) != d) stop("delta wrong dimension")
    out <- .C("astfam_cumulant", theta = as.double(theta), fam = as.integer(1),
        deriv = as.integer(deriv), delta = as.double(delta),
        zeroth = double(1), first = double(d),
        second = matrix(0, d, d), third = array(0, rep(d, 3)),
        PACKAGE = "aster2")
    fam.clear()
    result <- list(zeroth = out$zeroth)
    if (deriv >= 1) result$first <- out$first
    if (deriv >= 2) result$second <- out$second
    if (deriv >= 3) result$third <- out$third
    if (d == 1) result <- lapply(result, as.vector)
    return(result)
}

link <- function(xi, fam, deriv = 0, delta) {
    stopifnot(inherits(fam, "astfam"))
    stopifnot(is.atomic(xi))
    stopifnot(is.numeric(xi))
    stopifnot(all(is.finite(xi)))
    stopifnot(is.atomic(deriv))
    stopifnot(is.numeric(deriv))
    stopifnot(length(deriv) == 1)
    stopifnot(deriv == as.integer(deriv))
    stopifnot(deriv >= 0 && deriv <= 3)
    fam.clear()
    fam.set(fam)
    d <- fam.dimension(1)
    if (missing(delta)) delta <- rep(0, d)
    stopifnot(is.atomic(delta))
    stopifnot(is.numeric(delta))
    stopifnot(all(is.finite(delta)))
    if (length(xi) != d) stop("xi wrong dimension")
    if (length(delta) != d) stop("delta wrong dimension")
    out <- .C("astfam_link", xi = as.double(xi), fam = as.integer(1),
        deriv = as.integer(deriv), delta = as.double(delta),
        zeroth = double(d), first = matrix(0, d, d),
        PACKAGE = "aster2")
    fam.clear()
    result <- list(zeroth = out$zeroth)
    if (deriv >= 1) result$first <- out$first
    if (d == 1) result <- lapply(result, as.vector)
    return(result)
}

