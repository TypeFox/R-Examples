`fitted.timetrack` <-
    function(object, which = c("passive", "ordination"),
             model = NULL, choices = 1:2, ...)
{
    if(missing(which))
        which <- "passive"
    which <- match.arg(which)
    model <- if(is.null(model)) {
        if(is.null(object$ordination$CCA)) "CA" else "CCA"
    }
    if(isTRUE(all.equal(which, "passive"))) {
        fit <- fitted(unclass(object), ...)[, choices, drop = FALSE]
    } else {
        fit <- fitted(object$ordination, model = model,
                      ...)[, choices, drop = FALSE]
      }
    fit
}
