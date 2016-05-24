`scores.timetrack` <- function(x, which = c("ordination","passive"),
                               scaling = x$scaling, choices = 1:2,
                               display = "sites", ...) {
    which <- match.arg(which)
    ## only scores available for passive are the fitted (predicted) WA
    ## scores
    scrs <- if(which == "passive") {
        fitted(x, which = which, choices = choices, ...)
    } else {
        ## for the underlying ordination, take the usual scores
        scores(x$ord, ..., choices = choices, scaling = scaling,
               display = display)
    }
    scrs
}
