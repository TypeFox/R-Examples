
anova.asterOrReaster <- function(object, ...,
    tolerance = .Machine$double.eps ^ 0.75)
{
    dotargs <- list(...)
    if (length(dotargs) == 0)
        stop("need at least two models to compare")
    allargs <- c(list(object), dotargs)
    if (! all(sapply(allargs, function(x) inherits(x, "asterOrReaster"))))
        stop("some arguments not of class \"asterOrReaster\"")
    return(anovaAsterOrReasterList(allargs, tolerance))
}

anovaAsterOrReasterList <- function(objectlist,
    tolerance = .Machine$double.eps ^ 0.75)
{
    stopifnot(is.list(objectlist))
    if (length(objectlist) < 2)
        stop("must compare two or more models")

    if (! all(sapply(objectlist, function(x) inherits(x, "asterOrReaster"))))
        stop("some components not of class \"asterOrReaster\"")

    nmodels <- length(objectlist)

    stopifnot(is.numeric(tolerance))
    stopifnot(length(tolerance) == 1)
    stopifnot(tolerance > 0)

    # attempt to check that models are nested

    # FIXME: this is bad, but don't check offset vectors for now
    #     assume default offset vector was used

    get.fixed <- function(x) {
        if(inherits(x, "reaster")) {
            return(x$fixed)
        }
        if(inherits(x, "aster")) {
            bar <- x$modmat
            nparm <- dim(bar)[3]
            baz <- matrix(as.vector(bar), ncol = nparm)
            colnames(baz) <- dimnames(bar)[[3]]
            return(baz)
        }
        stop("model object not class \"aster\" or \"reaster\"")
    }
    fixed <- lapply(objectlist, get.fixed)

    get.random <- function(x) {
        if(inherits(x, "reaster"))
            return(x$random)
        if(inherits(x, "aster"))
            return(list())
        stop("model object not class \"aster\" or \"reaster\"")
    }
    random <- lapply(objectlist, get.random)

    nnode.fixed <- sapply(fixed, nrow)
    nnode.random <- lapply(fixed, function(x) sapply(x, nrow))
    nnode <- c(nnode.fixed, unlist(nnode.random))
    if (length(unique(nnode)) != 1)
        stop("model matrices do not all have the same row dimension")
    nnode <- unique(nnode)

    df.fixed <- sapply(fixed, ncol)
    df.random <- sapply(random, length)

    if (any(diff(df.fixed) < 0))
        stop("degrees of freedom for fixed effects not nondecreasing")
    if (any(diff(df.random) < 0))
        stop("degrees of freedom for random effects not nondecreasing")

    # check model matrix for fixed effects nesting

    qrList <- lapply(fixed, qr)
    for (i in 2:nmodels) {
        modmat1 <- fixed[[i - 1]]
        qr2 <- qrList[[i]]
        resid1 <- qr.resid(qr2, modmat1)
        norm.modmat1 <- apply(modmat1, 2, function(x) sqrt(sum(x^2)))
        norm.resid1 <- apply(resid1, 2, function(x) sqrt(sum(x^2)))
        if (any(norm.resid1 / norm.modmat1 > tolerance))
            stop("model matrices for fixed effects not nested\nmodel ",
                i - 1, " and model ", i)
    }

    # check model matrix for random effects nesting
    # here we hope simple check on column labels is o. k.

    collab.random <- lapply(random, function(x) lapply(x, colnames))
    foo.random <- lapply(collab.random, function(x) sapply(x, is.null))
    if (any(unlist(foo.random))) {
        warning("some model matrices for random effects have no colnames\nnot checking for nested models")
    } else {
        for (i in 2:nmodels) {
            c1 <- collab.random[[i - 1]]
            c2 <- collab.random[[i]]
            p1 <- sapply(c1, paste, collapse = " ")
            p2 <- sapply(c2, paste, collapse = " ")
            if (! all(p1 %in% p2))
                stop("model matrices for random effects not nested\nmodel ",
                i - 1, " and model ", i)
        }
    }

    resdf.fixed <- diff(df.fixed)
    resdf.random <- diff(df.random)

    if (any(resdf.random > 1))
        stop("don't know how to compare models differing by two or more variance components")

    resdev <- sapply(objectlist, function(x) x$deviance)

    get.formulae <- function(x) {
        if(inherits(x, "aster")) {
            qux <- x$formula
            if (is.null(qux))
                return("(no formulas)")
            return(as.character(deparse(qux)))
        }
        if(inherits(x, "reaster")) {
            qux <- x$formula
            if (is.null(qux))
                return("(no formulas)")
            qux <- unlist(qux)
            quacks <- sapply(qux, function(x) as.character(deparse(x)))
            return(paste(quacks, collapse = ", "))
        }
        stop("model object not class \"aster\" or \"reaster\"")
    }
    formulae <- sapply(objectlist, get.formulae)

    if (all(df.random == 0)) {
        # all aster models
        table <- data.frame(df.fixed, - resdev, c(NA, resdf.fixed),
            c(NA, diff(- resdev)))
        dimnames(table) <- list(1:nmodels,
            c("Model Df", "Model Dev", "Df", "Deviance"))
        title <- "Analysis of Deviance Table\n"
        topnote <- paste("Model ", format(1:nmodels), ": ",
            formulae, collapse = "\n", sep = "")
        table <- cbind(table, "P(>|Chi|)" = pchisq(table[ , "Deviance"],
            table[ , "Df"], lower.tail = FALSE))
        table <- structure(table, heading = c(title, topnote),
            class = c("anova", "data.frame"))
    } else {
        # some, perhaps all reaster models
        table <- data.frame(df.fixed, df.random, - resdev, c(NA, resdf.fixed),
            c(NA, resdf.random), c(NA, diff(- resdev)))
        dimnames(table) <- list(1:nmodels,
            c("Mod Df Fix", "Mod Df Rand", "Mod Dev",
            "Df Fix", "Df Rand", "Deviance"))
        title <- "Analysis of Deviance Table\n"
        topnote <- paste("Model ", format(1:nmodels), ": ",
            formulae, collapse = "\n", sep = "")
        d <- table[ , "Deviance"]
        nf <- table[ , "Df Fix"]
        nr <- table[ , "Df Rand"]
        p <- (1 / 2) * pchisq(d, nf, lower.tail = FALSE) +
            (1 / 2) * pchisq(d, nf + nr, lower.tail = FALSE)
        table <- cbind(table, "P-value" = p)
        table <- structure(table, heading = c(title, topnote),
            class = c("anova", "data.frame"))
    }

    return(table)
}

