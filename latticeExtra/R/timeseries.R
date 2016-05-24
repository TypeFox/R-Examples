
xyplot.stl <-
    function(x, data = NULL,
             outer = TRUE,
             layout = c(1, 4),
             strip = FALSE,
             strip.left = TRUE,
             as.table = TRUE,
             ylab = "",
             between = list(y = 0.5),
             panel =
             function(..., type) {
                 if (packet.number() == 4) type <- "h"
                 panel.xyplot(..., type = type)
             },
             ...)
{
    stopifnot(is.null(data))
    mstrip <- missing(strip.left)
    sers <- x$time.series
    ## ncomp <- ncol(sers)
    data <- rowSums(sers)
    X <- cbind(data, sers)
    colnames(X) <- c("data", colnames(sers))
    ans <-
        xyplot(X,
               outer = outer,
               layout = layout,
               strip = strip,
               strip.left = strip.left,
               as.table = as.table,
               ylab = ylab,
               between = between,
               panel = panel,
               ...,
               default.scales =
               list(x = list(axs = "i"),
                    y =
                    list(relation = "free",
                         tick.number = 3,
                         rot = 0)))
    if (mstrip)
    {
        mx <- min(rx <- abs(sapply(ans$y.limits, diff)))
        int <- cbind(-mx / rx, mx / rx)
        ans <-
            update(ans,
                   strip.left =
                   strip.custom(horizontal = FALSE,
                                strip.names = FALSE,
                                strip.levels = TRUE,
                                shingle.intervals = int))
    }
    ans$call <- sys.call(sys.parent()); ans$call[[1]] <- quote(xyplot)
    ans
}



