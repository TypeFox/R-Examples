`screeplot.pcr` <- function(x, restrict = NULL,
                            display = c("RMSE","avgBias","maxBias","R2"),
                            xlab = NULL, ylab = NULL, main = NULL, sub = NULL,
                            ...) {
    display <- match.arg(display)
    captions <- c("RMSE", "Average bias", "Maximum bias", "R squared")
    names(captions) <- c("RMSE", "avgBias", "maxBias", "R2")
    if (is.null(xlab))
        xlab <- "No. of components"
    if (is.null(ylab))
        ylab <- captions[display]
    if (is.null(main))
        main <- deparse(substitute(x))
    if (is.null(sub)) {
        cal <- x$call
        ##if (!is.na(m.f <- match("formula", names(cal)))) {
        ##    cal <- cal[c(1, m.f)]
        ##    names(cal)[2] <- ""
        ##}
        cc <- deparse(cal, 90)
        nc <- nchar(cc[1])
        abbr <- length(cc) > 1 || nc > 90
        sub <- if (abbr)
            paste(substr(cc[1], 1, min(90, nc)), "...")
        else cc[1]
    }
    dat <- performance(x)
    if(!is.null(restrict)) {
        comps <- min(restrict, x$ncomp)
    } else {
        comps <- x$ncomp
    }
    Scomps <- seq_len(comps)
    plot(Scomps, dat[Scomps, display], type = "n", ylab = ylab, xlab = xlab,
         main = main, sub = sub, ...)
    if(comps > 20) {
        lines(Scomps, dat[Scomps, display], type = "b", ...)
    } else {
        lines(Scomps, dat[Scomps, display], type = "b", ..., pch = NA)
        text(Scomps, dat[Scomps, display], labels = as.character(Scomps), cex = 0.8,
             ...)
    }
    invisible()
}
