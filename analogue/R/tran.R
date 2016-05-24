`tran` <- function(x, ...) {
    UseMethod("tran")
}

`tran.default` <- function(x, method, a = 1, b = 0, p = 2, base = exp(1),
                   na.rm = FALSE, na.value = 0, ...) {
    repMissing <- function(x, na.value) {
        x[is.na(x)] <- na.value
        return(x)
    }
    wasDF <- is.data.frame(x)
    dim.nams <- dimnames(x)
    x <- data.matrix(x)
    METHOD <- c("sqrt", "cubert", "rootroot", "log", "reciprocal", "freq",
                "center", "standardize", "range", "percent", "proportion",
                "pa","missing", "hellinger", "chi.square", "wisconsin",
                "pcent2prop", "prop2pcent", "logRatio", "power",
                "rowCentre", "colCentre", "rowCenter", "colCenter",
                "log1p", "expm1", "none")
    method <- match.arg(method, METHOD)
    ## account for non-British spelling
    american <- c("rowCenter", "colCenter")
    if (any(ind <-  american == method)) {
        method <- american[ind]
    }
    if(method %in% c("freq", "standardize","range","pa","hellinger",
                     "chi.square","wisconsin")) {
        if(isTRUE(all.equal(method, "wisconsin")))
            x <- wisconsin(x)
        else
            x <- decostand(x, method = method, na.rm = na.rm, ...)
        attr(x, "decostand") <- NULL
    } else {
        x <- switch(method,
                    sqrt = sqrt(x),
                    cubert = sign(x) * exp(log(abs(x)) / 3), #x^(1/3),
                    rootroot = sign(x) * exp(log(abs(x)) / 4), #x^(1/4),
                    log = {x <- sweep(x, 2, a, "*")
                           x <- sweep(x, 2, b, "+")
                           log(x, base = base)} ,
                    reciprocal = 1 / x,
                    center = scale(x, scale = FALSE, center = TRUE),
                    percent = sweep(x, 1, rowSums(x), "/") * 100,
                    proportion = sweep(x, 1, rowSums(x), "/"),
                    missing = repMissing(x, na.value),
                    pcent2prop = x / 100,
                    prop2pcent = x * 100,
                    logRatio = {x <- sweep(x, 2, a, "*")
                                x <- sweep(x, 2, b, "+")
                                x <- log(x, base = base)
                                x - rowMeans(x)},
                    power = x^p,
                    rowCentre = x - rowMeans(x),
                    colCentre = x - colMeans(x),
                    log1p = log1p(x),
                    expm1 = expm1(x),
                    none = x
                    )
    }
    if(wasDF)
        x <- as.data.frame(x)
    dimnames(x) <- dim.nams
    attr(x, "tran") <- method
    x
}

`tran.formula` <- function(formula, data = NULL,
                           subset = NULL,
                           na.action = na.pass, ...) {
    mf <- match.call()
    mf[[1]] <- as.name("model.frame")
    mt <- terms(formula, data = data, simplify = TRUE)
    mf[[2]] <- formula(mt, data = data)
    mf$na.action <- substitute(na.action)
    dots <- list(...)
    mf[names(dots)] <- NULL
    mf <- eval(mf, parent.frame())
    tran.default(mf, ...)
}
