`scaleFit` <-
function (fitted)
{
    ## m.call <- fitted$call
    ## mf <- model.frame(fitted$terms,data = eval(m.call$data))
    ## X.or <- model.matrix(fitted$terms, mf, contrasts = fitted$contrasts)
    ## if (inherits(fitted,"polr"))
    ##     X.or <- X.or[, -1]
    ## if (inherits(fitted,"BTm")) {
    ##     X.or <- fitted$x0
    ## }
    ## X.max.scaleFit <- apply(abs(X.or), 2, max)
    ## m.formula <- formula(fitted)
    ## offs <- model.offset(mf)
    ## if (is.null(offs))
    ##     offs <- rep(0, nrow(X.or))
    ## test1 <- ".the.scaled." %in% ls(.GlobalEnv)
    ## test2 <- ".the.offs." %in% ls(.GlobalEnv)
    ## assign(".the.scaled.", value = sweep(X.or, 2, X.max.scaleFit,
    ##     "/"), envir = .GlobalEnv)
    ## assign(".the.offs.", value = offs, envir = .GlobalEnv)
    ## browser()
    ## m.call$formula <- update.formula(m.formula, ~-1 + .the.scaled. +
    ##     offset(.the.offs.))
    ## m.call$offset <- NULL
    ## suppressWarnings(new.fit <- eval(m.call))
    ## new.fit$X.max.scaleFit <- X.max.scaleFit
    ## new.fit
    new.fit <- fitted
    new.fit$X.max.scaleFit <- 1
    new.fit
}
