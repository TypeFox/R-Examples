anova.logbin <- function(object, ..., test = NULL)
{
    dotargs <- list(...)
    named <- if (is.null(names(dotargs))) rep(FALSE, length(dotargs))
                else (names(dotargs) != "")
    if (any(named))
        warning("the following arguments to 'anova.logbin' are invalid and dropped: ",
            paste(deparse(dotargs[named]), collapse = ", "))
    dotargs <- dotargs[!named]
    is.logbin <- unlist(lapply(dotargs, function(x) inherits(x, "logbin")))
    dotargs <- dotargs[is.logbin]
    if (length(dotargs))
        return(anova.logbinlist(c(list(object), dotargs), test = test))
    else
        stop('anova.logbin does not support anova for a single model. Fit nested models manually and input to anova.logbin')
}