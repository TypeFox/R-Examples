`Within` <- function(type = c("free","series","grid","none"),
                     constant = FALSE, mirror = FALSE,
                     ncol = NULL, nrow = NULL)
{
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

    out <- list(type = type, constant = constant, mirror = mirror,
                ncol = ncol, nrow = nrow, call = .call)
    class(out) <- "Within"
    out
}
