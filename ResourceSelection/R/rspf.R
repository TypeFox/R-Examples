rspf <-
function(formula, data, m, B = 99, link = "logit", 
inits, method = "Nelder-Mead", control,
model = TRUE, x = FALSE, ...)
{
    link <- match.arg(link, c("logit","cloglog","probit","log"))
    if (link == "log")
        stop("link='log' is not allowed in rspf, use the rsf function instead")
    ## parsing formula
    if (missing(data))
        data <- parent.frame()
    mf <- match.call(expand.dots = FALSE)
    mm <- match(c("formula", "data"), names(mf), 0)
    mf <- mf[c(1, mm)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    Y <- model.response(mf, "numeric")
    ff <- formula
    ff[[2]] <- NULL
    mt <- terms(ff, data = data)
    X <- model.matrix(mt, mf)
    Xlevels <- .getXlevels(mt, mf)
    ## check variables
    if (length(Y) < 1) 
        stop("empty model")
    if (all(Y > 0)) 
        stop("invalid dependent variable, no zero value")
    if (!isTRUE(all.equal(as.vector(Y), as.integer(round(Y + 
        0.001))))) 
        stop("invalid dependent variable, non-integer values")
    Y <- as.integer(round(Y + 0.001))
    if (any(Y < 0)) 
        stop("invalid dependent variable, negative counts")
    if (any(!(Y %in% c(0, 1)))) 
        stop("invalid dependent variable, not in c(0, 1)")
    if (length(Y) != NROW(X)) 
        stop("invalid dependent variable, not a vector")
    if (identical(as.character(ff[[2]]), "1"))
        stop("invalid formula, no covariates")
    factonly <- all(unique(sapply(mf, .MFclass)[-1]) %in% c("ordered", "factor"))
    if (factonly)
        stop("provide at least 1 continuous covariate for RSPF")

    ## fitting
    out1 <- rsf.fit(X=X, Y=Y, m=m, link = link, B = B, 
        inits=inits, method = method, control=control, ...)

    ## return value assembly
    out1$call <- match.call()
    out2 <- list(formula=formula,
        terms=mt,
        levels=Xlevels,
        contrasts=attr(X, "contrasts"),
        model= if (model) mf else NULL,
        x= if (x) X else NULL)
    out <- c(out1, out2)
    ## defining object class
    class(out) <- c("rspf", "rsf")
    out
}

