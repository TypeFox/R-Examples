`how` <- function(within = Within(),
                  plots = Plots(),
                  blocks = NULL,
                  nperm = 199,
                  complete = FALSE,
                  maxperm = 9999,
                  minperm = 5040,
                  all.perms = NULL,
                  make = TRUE,
                  observed = FALSE) {

    blocks.name <- deparse(substitute(blocks))
    ## blocks should also be a factor - coerce
    if(!is.null(blocks))
        blocks <- as.factor(blocks)

    ## process the call to make it standalone
    .call <- match.call()
    if (length(.call) > 1L) {
        .ll <- as.list(.call[-1])
        args <- names(.call)[-1]
        ## evaluate arguments other than within and plots
        ## those handled in their respective functions
        for (i in args[!args %in% c("within","plots")]) {
            if(!is.null(.ll[[i]])) {
                .ll[[i]] <- eval(.ll[[i]], parent.frame())
            }
        }
    }

    out <- list(within = within, plots = plots, blocks = blocks,
                nperm = nperm, complete = complete,
                maxperm = maxperm, minperm = minperm,
                all.perms = all.perms, make = make,
                observed = observed,
                blocks.name = blocks.name)

    ## process within and plots separately
    if (length(.call) > 1L && "within" %in% args) {
        .ll[["within"]] <- getCall(within)
    }
    if (length(.call) > 1L && "plots" %in% args) {
        .ll[["plots"]] <- getCall(plots)
    }

    ## finsh off
    if (length(.call) > 1L) {
        .ll <- c(as.list(.call[[1]]), .ll)
        names(.ll) <- names(.call)
        .call <- as.call(.ll)
    }

    out$call <- .call

    class(out) <- "how"
    out
}
