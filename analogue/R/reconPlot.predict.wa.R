reconPlot.predict.wa <-
    function(x, depths, use.labels = FALSE,
             display.error = c("none","bars","lines"),
             sample.specific = TRUE, ...) {
        if (missing(display.error))
            display.error <- "none"
        display.error <- match.arg(display.error)
        if (missing(depths)) {
            if (use.labels)
                depths <- as.numeric(names(x$pred$pred))
            else
                stop("If \"use.labels = FALSE\", then \"depths\" must be provided.")
        }
        preds <- x$pred$pred
        if(sample.specific && !is.null(x$pred$rmsep))
            errors <- x$pred$rmsep
        else
            errors <- rep(x$performance$rmsep, length(x$pred$pred))
        reconPlot.default(preds, depths, errors,
                          display.error = display.error, ...)
    }
