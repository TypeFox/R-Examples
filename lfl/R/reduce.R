reduce <- function(x,
                    rules,
                    ratio,
                    tnorm=c("goedel", "goguen", "lukasiewicz"),
                    tconorm=c("goedel", "goguen", "lukasiewicz"),
                    numThreads=1) {
    if (is.vector(x)) {
        x <- matrix(x, nrow=1, dimnames=list(NULL, names(x)))
    }
    if (!is.matrix(x) || !is.numeric(x)) {
        stop("'x' must be a numeric vector or matrix")
    }
    if (min(x) < 0 || max(x) > 1) {
        stop("Values of 'x' must be truth values in the interval [0,1]")
    }

    origRules <- rules
    if (is.farules(origRules)) {
        rules <- origRules$rules
    } else if (is.vector(rules) && is.character(rules)) {
        rules <- list(rules)
    }
    if (!is.list(rules)) {
        stop("'rules' must be a list of rules")
    }

    if (!is.numeric(numThreads) || length(numThreads) > 1 || numThreads < 1) {
        stop("'numThreads' must be positive integer number")
    }

    tnorm <- match.arg(tnorm)
    if (tnorm == 'goedel') {
        tnorm <- 'minimum'
    } else if (tnorm == 'goguen') {
        tnorm <- 'product'
    }

    tconorm <- match.arg(tconorm)
    if (tconorm == 'goedel') {
        tconorm <- 'maximum'
    } else if (tconorm == 'goguen') {
        tconorm <- 'product'
    }

    lhsSupport <- NA
    if (is.farules(origRules) && ('lhsSupport' %in% colnames(origRules$statistics))) {
        lhsSupport <- origRules$statistics[, 'lhsSupport']
    } else {
        lhsSupport <- rep(1, length(rules))
    }

    lvl <- colnames(x)
    rb <- lapply(rules, function(x) { as.integer(factor(x, levels=lvl)) - 1 })
    config <- list(data=x,
                   rules=rb,
                   lhsSupport=lhsSupport,
                   ratio=ratio,
                   tnorm=tnorm,
                   tconorm=tconorm,
                   numThreads=numThreads)
    result <- .Call("reduce", config, PACKAGE="lfl")
    result <- result + 1

    if (is.farules(origRules)) {
        r <- origRules$rules[result]
        s <- origRules$statistics[result, , drop=FALSE]
        return(farules(rules=r, statistics=s))
    } else {
        return(rules[result])
    }
}
