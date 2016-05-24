model.frame.intReg <- function (formula, ...) {
   ## 'formula' is supposed to be an  'intReg' object
    dots <- list(...)
    nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 0 )]
    if ( length(nargs) > 0 || is.null(formula$model)) {
        fcall <- formula$call
        fcall$method <- "model.frame"
        fcall[[1]] <- as.name("intReg")
        fcall[names(nargs)] <- nargs
        env <- environment(formula$terms)
        if (is.null(env)) {
            env <- parent.frame()
        }
        result <- eval(fcall, env, parent.frame())
    } else {
      result <- formula$model
    }
    return( result )
}
