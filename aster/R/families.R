
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

fam.truncated.poisson <- function(truncation) {
    stopifnot(is.numeric(truncation))
    stopifnot(length(truncation) == 1)
    stopifnot(truncation == round(truncation))
    stopifnot(truncation >= 0)
    result <- list(name = "truncated.poisson", truncation = truncation)
    class(result) <- "astfam"
    return(result)
}

fam.negative.binomial <- function(size) {
    stopifnot(is.numeric(size))
    stopifnot(length(size) == 1)
    stopifnot(size >= 0)
    result <- list(name = "negative.binomial", size = size)
    class(result) <- "astfam"
    return(result)
}

fam.truncated.negative.binomial <- function(size, truncation) {
    stopifnot(is.numeric(size))
    stopifnot(length(size) == 1)
    stopifnot(size >= 0)
    stopifnot(is.numeric(truncation))
    stopifnot(length(truncation) == 1)
    stopifnot(truncation == round(truncation))
    stopifnot(truncation >= 0)
    result <- list(name = "truncated.negative.binomial", size = size,
        truncation = truncation)
    class(result) <- "astfam"
    return(result)
}

fam.normal.location <- function(sd) {
    stopifnot(is.numeric(sd))
    stopifnot(length(sd) == 1)
    stopifnot(sd > 0)
    result <- list(name = "normal.location", sd = sd)
    class(result) <- "astfam"
    return(result)
}

fam.default <- function() {
    list(fam.bernoulli(), fam.poisson(), fam.truncated.poisson(truncation = 0),
       fam.truncated.poisson(truncation = 2))
}

clearfam <- function()
    invisible(.C("aster_clear_families", PACKAGE = "aster"))

getsupfambyname <- function(fam) {
    stopifnot(is.character(fam))
    stopifnot(length(fam) == 1)
    foo <- .C("aster_byname_superfamily", name = fam, nhyper = integer(1),
        hypername = character(2), PACKAGE = "aster")
    foo$name <- fam
    return(foo)
}

setfam <- function(famlist) {
    stopifnot(is.list(famlist))
    stopifnot(all(sapply(famlist, inherits, what = "astfam")))
    clearfam()
    for (i in seq(along = famlist)) {
        fam <- famlist[[i]]
        famname <- fam$name
        foo <- getsupfambyname(famname)
        hyper <- double(2)
        if (foo$nhyper >= 1) {
            qux <- foo$hypername[1]
            bar <- fam[[qux]]
            if (is.null(bar))
                stop("family \"", famname, "\" needs hyperparameter \"",
                    foo$hypername[1], "\"")
            hyper[1] <- bar
        }
        if (foo$nhyper >= 2) {
            qux <- foo$hypername[2]
            bar <- fam[[qux]]
            if (is.null(bar))
                stop("family \"", famname, "\" needs hyperparameter \"",
                    foo$hypername[2], "\"")
            hyper[2] <- bar
        }
       .C("aster_add_family", name = famname, hyper = as.double(hyper),
           nhyper = as.integer(foo$nhyper), PACKAGE = "aster")
    }
}

getfam <- function() {
    result <- list()
    ifam <- 1
    repeat {
        foo <- .C("aster_get_family", idx = as.integer(ifam),
            name = character(1), hyper = double(2), nhyper = integer(1),
            hypername = character(2), origin = double(1), PACKAGE = "aster")
        if (foo$name == "")
            break;
        sally <- list(name = foo$name)
        if (foo$nhyper >= 1) {
            bar <- foo$hypername[1]
            baz <- foo$hyper[1]
            sally[[bar]] <- baz
        }
        if (foo$nhyper >= 2) {
            bar <- foo$hypername[2]
            baz <- foo$hyper[2]
            sally[[bar]] <- baz
        }
        sally$origin <- foo$origin
        result[[ifam]] <- sally
        ifam <- ifam + 1
    }
    return(result)
}

getsupfam <- function() {
    result <- list()
    ifam <- 1
    repeat {
        foo <- .C("aster_get_superfamily", idx = as.integer(ifam),
            name = character(1), nhyper = integer(1),
            hypername = character(2), PACKAGE = "aster")
        if (foo$name == "")
            break;
        sally <- list(name = foo$name)
        if (foo$nhyper >= 1) {
            bar <- foo$hypername[1]
            sally[[bar]] <- NA
        }
        if (foo$nhyper >= 2) {
            bar <- foo$hypername[2]
            sally[[bar]] <- NA
        }
        result[[ifam]] <- sally
        ifam <- ifam + 1
    }
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
    print(foo)
    return(invisible(foo))
}

