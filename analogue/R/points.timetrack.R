`points.timetrack` <- function(x, choices = 1:2,
                               which = c("passive", "ordination"),
                               display = c("wa","lc"),
                               order,
                               ...) {
    display <- match.arg(display)
    which <- match.arg(which)

    ## Select the coordinates for the relevant type of sample
    if (isTRUE(all.equal(which, "ordination"))) {
        scrs <- scores(x$ord, choices = choices, scaling = x$scaling,
                       display = display, ...)
    } else {
        scrs <- fitted(x, type = "passive", choices = choices)
        if(!missing(order)) {
            if(length(order) != NROW(scrs))
                stop("'length(order)' not equal to number of passive samples.")
            scrs[order, ]
        }
    }

    points(scrs, ...)
}
