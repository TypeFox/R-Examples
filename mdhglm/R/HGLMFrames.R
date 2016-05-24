HGLMFrames <-
function (mc, formula, contrasts, vnms = character(0)) 
{
    mf <- mc
    m <- match(c("DataMain", "weights", "na.action", "offset"), 
        names(mf), 0)
    mf <- mf[c(1, m)]
    frame.form <- subbars(formula)
    if (length(vnms) > 0) 
        frame.form[[3]] <- substitute(foo + bar, list(foo = parse(text = paste(vnms, 
            collapse = " + "))[[1]], bar = frame.form[[3]]))
    fixed.form <- nobars(formula)
    if (inherits(fixed.form, "name")) 
        fixed.form <- substitute(foo ~ 1, list(foo = fixed.form))
    environment(fixed.form) <- environment(frame.form) <- environment(formula)
    mf$formula <- frame.form
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    names(mf)[2] <- "data"
    fe <- mf
    mf <- eval(mf, parent.frame(2))
    fe$formula <- fixed.form
    fe <- eval(fe, parent.frame(2))
    fe
    Y <- model.response(mf, "any")
    if (length(dim(Y)) == 1) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if (!is.null(nm)) 
            names(Y) <- nm
    }
    mt <- attr(fe, "terms")
    X <- if (!is.empty.model(mt)) 
        model.matrix(mt, mf, contrasts)
    else matrix(, NROW(Y), 0)
    storage.mode(X) <- "double"
    fixef <- numeric(ncol(X))
    names(fixef) <- colnames(X)
    dimnames(X) <- NULL
    wts <- model.weights(mf)
    if (is.null(wts)) 
        wts <- numeric(0)
    off <- model.offset(mf)
    if (is.null(off)) 
        off <- numeric(0)
    if (any(wts <= 0)) 
        stop(gettextf("negative weights or weights of zero are not allowed"))
    if (length(off) && length(off) != NROW(Y)) 
        stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
            length(off), NROW(Y)))
    attr(mf, "terms") <- mt
    list(Y = Y, X = X, wts = as.double(wts), off = as.double(off), 
        mf = mf, fixef = fixef)
}
