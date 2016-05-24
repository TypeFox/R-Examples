`plot3d.prcurve` <- function(x, choices = 1:3, display = "sites",
                             scaling = 0,
                             lcol = "darkorange", lwd = 2,
                             decorate = TRUE,
                             xlab = NULL, ylab = NULL, zlab = NULL,
                             main = NULL, ...) {
    ## if (!require("rgl")) {
    ##     stop("Package 'rgl' required but not installed.\nRun\n\t'install.packages(\"rgl\")\nto use 'Plot3d()'.")
    ## }
    ## check choices
    if(!missing(choices)) {
        lc <- length(choices)
        if(lc > 3L) {
            warning("More than 3 axes specified in 'choices', using the first 3.")
            choices <- rep(choices, length.out = 3)
        } else if(lc < 3L) {
            warning("Fewer than 3 axes specified; reverting to defaults ('1:3')")
            choices <- 1:3
        }
    }
    ## do ordination
    ord <- x$ordination ## this is now stored
    ## process labels
    if(missing(xlab) || is.null(xlab))
        xlab <- paste0("PC", choices[1])
    if(missing(ylab) || is.null(ylab))
        ylab <- paste0("PC", choices[2])
    if(missing(zlab) || is.null(zlab))
        zlab <- paste0("PC", choices[3])
    ordirgl(ord, scaling = scaling, choices = choices, display = display,
            ...)
    pred <- predict(ord, x[["s"]], type = "wa",
                    scaling = scaling)[x[["tag"]], ]
    lines3d(pred[,1], pred[,2], pred[,3], col = lcol, lwd = lwd, ...)
    if(decorate) {
        decorate3d(xlab = xlab, ylab = ylab, zlab = zlab, ...)
    }
    invisible()
}
