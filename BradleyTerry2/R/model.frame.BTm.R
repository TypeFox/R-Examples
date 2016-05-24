model.frame.BTm <- function (formula, ...)
{
    dots <- list(...)
    nargs <- dots[match(c("outcome", "player1", "player2",
                          "separate.ability", "refcat", "data", "weights",
                          "subset", "offset", "contrasts"), names(dots), 0L)]
    mfArgs <- dots[match(c("na.action", "start", "etastart", "mustart"),
                         names(dots), 0L)]
    if (length(nargs) || is.null(formula$model)) {
        fcall <- formula$call[-1]
        fcall[names(nargs)] <- nargs
        env <- environment(formula$terms)
        if (is.null(env))
            env <- parent.frame()
        setup <- do.call(BTm.setup, fcall, envir = env)
        mf <- data.frame(X = setup$X[,1])
        mf$X <- setup$X
        mf$Y <- setup$Y
        mf <- as.call(c(model.frame, mfArgs,
                         list(formula = Y ~ X - 1, data = mf,
                              offset = setup$offset,
                              subset = setup$subset,
                              weights = setup$weights)))
        eval(mf, parent.frame())
    }
    else formula$model
}
