`Plots` <- function(strata = NULL, type = c("none","free","series","grid"),
                    mirror = FALSE, ncol = NULL, nrow = NULL) {

    plots.name <- deparse(substitute(strata))
    ## strata should also be a factor - coerce
    if(!is.null(strata))
        strata <- as.factor(strata)

    type <- match.arg(type)

    ## process the call to make it standalone
    .call <- match.call()
    if (length(.call) > 1L) {
        .ll <- as.list(.call[-1])
        for (i in seq_along(.ll))
            .ll[[i]] <- eval(.ll[[i]], parent.frame())
        .ll <- c(as.list(.call[[1]]), .ll)
        names(.ll) <- names(.call)
        .call <- as.call(.ll)
    }

    out <- list(strata = strata, type = type, mirror = mirror,
                ncol = ncol, nrow = nrow,
                plots.name = plots.name, call = .call)
    class(out) <- "Plots"
    out
}
