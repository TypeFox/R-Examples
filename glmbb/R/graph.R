findCliques <- function(mt) {

    stopifnot(inherits(mt, "terms"))

    mf <- attr(mt, "factors")
    mr <- attr(mt, "response")
    if (! is.matrix(mf)) return(list())
    if (mr != 0)
        mf <- mf[- mr, , drop = FALSE]
    g <- matrix(0, nrow(mf), nrow(mf))
    for (i in 1:nrow(mf))
        for (j in 1:nrow(mf))
            if (i != j)
                g[i, j] <- any(mf[i, ] * mf[j, ] >= 1)

    e <- new.env(hash = TRUE, parent = emptyenv())

    BronKerbosch <- function(R, P, X) {
        stopifnot(R %in% 1:nrow(mf))
        stopifnot(P %in% 1:nrow(mf))
        stopifnot(X %in% 1:nrow(mf))
        if (length(P) == 0 && length(X) == 0) {
           varname <- paste(c("foo", R), collapse = ".")
           assign(varname, R, envir = e)
        }
        for (v in P) {
            nv <- which(g[v, ] == 1)
            BronKerbosch(union(R, v), intersect(P, nv), intersect(X, nv))
            P <- setdiff(P, v)
            X <- union(X, v)
        }
    }

    foo <- BronKerbosch(integer(0), 1:nrow(mf), integer(0))
    bar <- as.list(e)
    names(bar) <- NULL
    bar
}

isGraphical <- function(formula) {

    stopifnot(inherits(formula, "formula"))

    mt <- terms(formula)
    mf <- attr(mt, "factors")
    if (! is.matrix(mf))
        return(TRUE)
    mr <- attr(mt, "response")
    if (mr != 0)
        mf <- mf[- mr, , drop = FALSE]

    if (! is.hierarchical(mt))
        stop("model is not hierarchical")

    bar <- findCliques(mt)

    ok <- TRUE
    for (i in bar) {
        qux <- 1:nrow(mf) %in% i
        quux <- apply(mf, 2, function(x) all((x >= 1) == qux))
        ok <- ok & any(quux)
    }
    return(ok)
}

asGraphical <- function(formula) {

    stopifnot(inherits(formula, "formula"))

    mt <- terms(formula)
    mf <- attr(mt, "factors")
    mr <- attr(mt, "response")
    mi <- attr(mt, "intercept")
    if (! is.matrix(mf))
        return(formula)

    if (mr != 0) {
        response.name <- rownames(mf)[mr]
        mf <- mf[- mr, , drop = FALSE]
    } else {
        response.name <- ""
    }
    covariates <- rownames(mf)

    bar <- findCliques(mt)
    bar <- lapply(bar, function(i) covariates[i])
    bar <- lapply(bar, function(x) paste(x, collapse = "*"))
    bar <- unlist(bar)
    bar <- paste(bar, collapse = " + ")
    bar <- paste(response.name, "~", bar)
    return(as.formula(bar))
}

isHierarchical <- function(formula) {
    stopifnot(inherits(formula, "formula"))
    is.hierarchical(terms(formula))
}

asHierarchical <- function(formula) {
    stopifnot(inherits(formula, "formula"))
    mt <- terms(formula)
    ml <- attr(mt, "term.labels")
    mi <- attr(mt, "intercept")
    mr <- attr(mt, "response")
    mv <- attr(mt, "variables")
    foo <- gsub(":", "*", ml)
    if (length(foo) == 0) {
        foo <- as.character(mi)
    } else {
        foo <- paste(foo, collapse = " + ")
        if (mi == 0)
            foo <- paste("0 +", foo)
    }
    foo <- paste("~", foo)
    if (mr != 0) {
        bar <- as.list(mv)[[mr + 1]]
        foo <- paste(bar, foo)
    }
    as.formula(foo)
}

