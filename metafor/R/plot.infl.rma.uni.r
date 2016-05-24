plot.infl.rma.uni <-
function (x, plotinf = TRUE, plotdfb = FALSE, dfbnew = FALSE, 
    logcov = TRUE, layout, slab.style = 1, las = 0, pch = 21, 
    bg = "black", bg.infl = "red", col.na = "lightgray", ...) 
{
    if (class(x) != "infl.rma.uni") 
        stop("Argument 'x' must be an object of class \"infl.rma.uni\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    any.na <- is.na(cbind(x$inf, x$dfb))
    if (any(any.na)) {
        if (na.act == "na.fail") 
            stop("Missing values in results.")
    }
    if (is.logical(plotinf)) {
        if (plotinf) {
            which.inf <- 1:8
        }
    }
    else {
        which.inf <- plotinf
        which.inf <- which.inf[(which.inf >= 1) & (which.inf <= 
            8)]
        which.inf <- unique(round(which.inf))
        if (length(which.inf) == 0L) 
            stop("Incorrect specification of 'plotinf' argument.")
        plotinf <- TRUE
    }
    if (is.logical(plotdfb)) {
        if (plotdfb) {
            which.dfb <- seq_len(x$p)
        }
    }
    else {
        which.dfb <- plotdfb
        which.dfb <- which.dfb[(which.dfb >= 1) & (which.dfb <= 
            x$p)]
        which.dfb <- unique(round(which.dfb))
        if (length(which.dfb) == 0L) 
            stop("Incorrect specification of 'plotdfb' argument.")
        plotdfb <- TRUE
    }
    if (!plotinf & !plotdfb) 
        stop("At least one of the arguments 'plotinf' or 'plotdfb' argument must be TRUE.")
    if (!plotinf & dfbnew) 
        dfbnew <- FALSE
    par.mar <- par("mar")
    par.mar.adj <- par.mar - c(2, 2, 2, 1)
    par.mar.adj[par.mar.adj < 1] <- 1
    par(mar = par.mar.adj)
    on.exit(par(mar = par.mar))
    lplot <- function(..., minlength, strict) {
        plot(...)
    }
    lpoints <- function(..., minlength, strict) {
        points(...)
    }
    llines <- function(..., minlength, strict) {
        lines(...)
    }
    laxis <- function(..., minlength, strict) {
        axis(...)
    }
    labline <- function(..., minlength, strict) {
        abline(...)
    }
    ids <- switch(slab.style, `1` = x$ids, `2` = rownames(x$inf), 
        `3` = abbreviate(rownames(x$inf), ...))
    if (plotinf) {
        par.mfrow <- par("mfrow")
        on.exit(par(mfrow = par.mfrow), add = TRUE)
        if (missing(layout)) {
            if (length(which.inf) == 2) 
                par(mfrow = c(2, 1))
            if (length(which.inf) == 3) 
                par(mfrow = c(3, 1))
            if (length(which.inf) == 4) 
                par(mfrow = c(2, 2))
            if (length(which.inf) == 5) 
                par(mfrow = c(5, 1))
            if (length(which.inf) == 6) 
                par(mfrow = c(3, 2))
            if (length(which.inf) == 7) 
                par(mfrow = c(7, 1))
            if (length(which.inf) == 8) 
                par(mfrow = c(4, 2))
        }
        else {
            layout <- layout[(layout >= 1)]
            layout <- round(layout)
            if (length(layout) != 2L) 
                stop("Incorrect specification of 'layout' argument.")
            par(mfrow = layout)
        }
        for (i in seq_len(length(which.inf))) {
            if (which.inf[i] == 1) {
                zi <- x$inf$rstudent
                not.na <- !is.na(zi)
                if (na.act == "na.omit") {
                  zi <- zi[not.na]
                  len.ids <- length(x$ids) - sum(!not.na)
                  ids.infl <- x$is.infl[not.na]
                  lab.ids <- ids[not.na]
                }
                if (na.act == "na.exclude" || na.act == "na.pass") {
                  len.ids <- length(x$ids)
                  ids.infl <- x$is.infl
                  lab.ids <- ids
                }
                zi.min <- min(zi, -2, na.rm = TRUE)
                zi.max <- max(zi, 2, na.rm = TRUE)
                lplot(NA, NA, xlim = c(1, len.ids), ylim = c(zi.min, 
                  zi.max), xaxt = "n", main = "rstudent", xlab = "", 
                  ylab = "", las = las, ...)
                laxis(side = 1, at = seq_len(len.ids), labels = lab.ids, 
                  xlab = "", las = las, ...)
                labline(h = 0, lty = "dashed", ...)
                labline(h = c(qnorm(0.025), qnorm(0.975)), lty = "dotted", 
                  ...)
                if (na.act == "na.exclude" || na.act == "na.pass") 
                  llines(seq_len(len.ids)[not.na], zi[not.na], 
                    col = col.na, ...)
                llines(seq_len(len.ids), zi, ...)
                lpoints(seq_len(len.ids), zi, pch = pch, bg = bg, 
                  ...)
                lpoints(seq_len(len.ids)[ids.infl], zi[ids.infl], 
                  bg = bg.infl, pch = pch, ...)
            }
            if (which.inf[i] == 2) {
                zi <- x$inf$dffits
                not.na <- !is.na(zi)
                if (na.act == "na.omit") {
                  zi <- zi[not.na]
                  len.ids <- length(x$ids) - sum(!not.na)
                  ids.infl <- x$is.infl[not.na]
                  lab.ids <- ids[not.na]
                }
                if (na.act == "na.exclude" || na.act == "na.pass") {
                  len.ids <- length(x$ids)
                  ids.infl <- x$is.infl
                  lab.ids <- ids
                }
                zi.min <- min(zi, na.rm = TRUE)
                zi.max <- max(zi, na.rm = TRUE)
                lplot(NA, NA, xlim = c(1, len.ids), ylim = c(zi.min, 
                  zi.max), xaxt = "n", main = "dffits", xlab = "", 
                  ylab = "", las = las, ...)
                laxis(side = 1, at = seq_len(len.ids), labels = lab.ids, 
                  xlab = "", las = las, ...)
                labline(h = 0, lty = "dashed", ...)
                labline(h = 3 * sqrt(x$p/(x$k - x$p)), lty = "dotted", 
                  ...)
                labline(h = -3 * sqrt(x$p/(x$k - x$p)), lty = "dotted", 
                  ...)
                if (na.act == "na.exclude" || na.act == "na.pass") 
                  llines(seq_len(len.ids)[not.na], zi[not.na], 
                    col = col.na, ...)
                llines(seq_len(len.ids), zi, ...)
                lpoints(seq_len(len.ids), zi, pch = pch, bg = bg, 
                  ...)
                lpoints(seq_len(len.ids)[ids.infl], zi[ids.infl], 
                  bg = bg.infl, pch = pch, ...)
            }
            if (which.inf[i] == 3) {
                zi <- x$inf$cook.d
                not.na <- !is.na(zi)
                if (na.act == "na.omit") {
                  zi <- zi[not.na]
                  len.ids <- length(x$ids) - sum(!not.na)
                  ids.infl <- x$is.infl[not.na]
                  lab.ids <- ids[not.na]
                }
                if (na.act == "na.exclude" || na.act == "na.pass") {
                  len.ids <- length(x$ids)
                  ids.infl <- x$is.infl
                  lab.ids <- ids
                }
                zi.min <- min(zi, na.rm = TRUE)
                zi.max <- max(zi, na.rm = TRUE)
                lplot(NA, NA, xlim = c(1, len.ids), ylim = c(zi.min, 
                  zi.max), xaxt = "n", main = "cook.d", xlab = "", 
                  ylab = "", las = las, ...)
                laxis(side = 1, at = seq_len(len.ids), labels = lab.ids, 
                  xlab = "", las = las, ...)
                labline(h = qchisq(0.5, df = x$p), lty = "dotted", 
                  ...)
                if (na.act == "na.exclude" || na.act == "na.pass") 
                  llines(seq_len(len.ids)[not.na], zi[not.na], 
                    col = col.na, ...)
                llines(seq_len(len.ids), zi, ...)
                lpoints(seq_len(len.ids), zi, pch = pch, bg = bg, 
                  ...)
                lpoints(seq_len(len.ids)[ids.infl], zi[ids.infl], 
                  bg = bg.infl, pch = pch, ...)
            }
            if (which.inf[i] == 4) {
                zi <- x$inf$cov.r
                not.na <- !is.na(zi)
                if (na.act == "na.omit") {
                  zi <- zi[not.na]
                  len.ids <- length(x$ids) - sum(!not.na)
                  ids.infl <- x$is.infl[not.na]
                  lab.ids <- ids[not.na]
                }
                if (na.act == "na.exclude" || na.act == "na.pass") {
                  len.ids <- length(x$ids)
                  ids.infl <- x$is.infl
                  lab.ids <- ids
                }
                zi.min <- min(zi, na.rm = TRUE)
                zi.max <- max(zi, na.rm = TRUE)
                if (logcov) {
                  lplot(NA, NA, xlim = c(1, len.ids), ylim = c(zi.min, 
                    zi.max), xaxt = "n", main = "cov.r", xlab = "", 
                    ylab = "", las = las, log = "y", ...)
                }
                else {
                  lplot(NA, NA, xlim = c(1, len.ids), ylim = c(zi.min, 
                    zi.max), xaxt = "n", main = "cov.r", xlab = "", 
                    ylab = "", las = las, ...)
                }
                laxis(side = 1, at = seq_len(len.ids), labels = lab.ids, 
                  xlab = "", las = las, ...)
                labline(h = 1, lty = "dashed", ...)
                if (na.act == "na.exclude" || na.act == "na.pass") 
                  llines(seq_len(len.ids)[not.na], zi[not.na], 
                    col = col.na, ...)
                llines(seq_len(len.ids), zi, ...)
                lpoints(seq_len(len.ids), zi, pch = pch, bg = bg, 
                  ...)
                lpoints(seq_len(len.ids)[ids.infl], zi[ids.infl], 
                  bg = bg.infl, pch = pch, ...)
            }
            if (which.inf[i] == 5) {
                zi <- x$inf$tau2.del
                not.na <- !is.na(zi)
                if (na.act == "na.omit") {
                  zi <- zi[not.na]
                  len.ids <- length(x$ids) - sum(!not.na)
                  ids.infl <- x$is.infl[not.na]
                  lab.ids <- ids[not.na]
                }
                if (na.act == "na.exclude" || na.act == "na.pass") {
                  len.ids <- length(x$ids)
                  ids.infl <- x$is.infl
                  lab.ids <- ids
                }
                zi.min <- min(zi, na.rm = TRUE)
                zi.max <- max(zi, na.rm = TRUE)
                lplot(NA, NA, xlim = c(1, len.ids), ylim = c(zi.min, 
                  zi.max), xaxt = "n", main = "tau2.del", xlab = "", 
                  ylab = "", las = las, ...)
                laxis(side = 1, at = seq_len(len.ids), labels = lab.ids, 
                  xlab = "", las = las, ...)
                labline(h = x$tau2, lty = "dashed", ...)
                if (na.act == "na.exclude" || na.act == "na.pass") 
                  llines(seq_len(len.ids)[not.na], zi[not.na], 
                    col = col.na, ...)
                llines(seq_len(len.ids), zi, ...)
                lpoints(seq_len(len.ids), zi, pch = pch, bg = bg, 
                  ...)
                lpoints(seq_len(len.ids)[ids.infl], zi[ids.infl], 
                  bg = bg.infl, pch = pch, ...)
            }
            if (which.inf[i] == 6) {
                zi <- x$inf$QE.del
                not.na <- !is.na(zi)
                if (na.act == "na.omit") {
                  zi <- zi[not.na]
                  len.ids <- length(x$ids) - sum(!not.na)
                  ids.infl <- x$is.infl[not.na]
                  lab.ids <- ids[not.na]
                }
                if (na.act == "na.exclude" || na.act == "na.pass") {
                  len.ids <- length(x$ids)
                  ids.infl <- x$is.infl
                  lab.ids <- ids
                }
                zi.min <- min(zi, na.rm = TRUE)
                zi.max <- max(zi, na.rm = TRUE)
                lplot(NA, NA, xlim = c(1, len.ids), ylim = c(zi.min, 
                  zi.max), xaxt = "n", main = "QE.del", xlab = "", 
                  ylab = "", las = las, ...)
                laxis(side = 1, at = seq_len(len.ids), labels = lab.ids, 
                  xlab = "", las = las, ...)
                labline(h = x$QE, lty = "dashed", ...)
                if (na.act == "na.exclude" || na.act == "na.pass") 
                  llines(seq_len(len.ids)[not.na], zi[not.na], 
                    col = col.na, ...)
                llines(seq_len(len.ids), zi, ...)
                lpoints(seq_len(len.ids), zi, pch = pch, bg = bg, 
                  ...)
                lpoints(seq_len(len.ids)[ids.infl], zi[ids.infl], 
                  bg = bg.infl, pch = pch, ...)
            }
            if (which.inf[i] == 7) {
                zi <- x$inf$hat
                not.na <- !is.na(zi)
                if (na.act == "na.omit") {
                  zi <- zi[not.na]
                  len.ids <- length(x$ids) - sum(!not.na)
                  ids.infl <- x$is.infl[not.na]
                  lab.ids <- ids[not.na]
                }
                if (na.act == "na.exclude" || na.act == "na.pass") {
                  len.ids <- length(x$ids)
                  ids.infl <- x$is.infl
                  lab.ids <- ids
                }
                zi.min <- 0
                zi.max <- max(zi, na.rm = TRUE)
                lplot(NA, NA, xlim = c(1, len.ids), ylim = c(zi.min, 
                  zi.max), xaxt = "n", main = "hat", xlab = "", 
                  ylab = "", las = las, ...)
                laxis(side = 1, at = seq_len(len.ids), labels = lab.ids, 
                  xlab = "", las = las, ...)
                labline(h = x$p/x$k, lty = "dashed", ...)
                labline(h = 3 * x$p/x$k, lty = "dotted", ...)
                if (na.act == "na.exclude" || na.act == "na.pass") 
                  llines(seq_len(len.ids)[not.na], zi[not.na], 
                    col = col.na, ...)
                llines(seq_len(len.ids), zi, ...)
                lpoints(seq_len(len.ids), zi, pch = pch, bg = bg, 
                  ...)
                lpoints(seq_len(len.ids)[ids.infl], zi[ids.infl], 
                  bg = bg.infl, pch = pch, ...)
            }
            if (which.inf[i] == 8) {
                zi <- x$inf$weight
                not.na <- !is.na(zi)
                if (na.act == "na.omit") {
                  zi <- zi[not.na]
                  len.ids <- length(x$ids) - sum(!not.na)
                  ids.infl <- x$is.infl[not.na]
                  lab.ids <- ids[not.na]
                }
                if (na.act == "na.exclude" || na.act == "na.pass") {
                  len.ids <- length(x$ids)
                  ids.infl <- x$is.infl
                  lab.ids <- ids
                }
                zi.min <- 0
                zi.max <- max(zi, na.rm = TRUE)
                lplot(NA, NA, xlim = c(1, len.ids), ylim = c(zi.min, 
                  zi.max), xaxt = "n", main = "weight", xlab = "", 
                  ylab = "", las = las, ...)
                laxis(side = 1, at = seq_len(len.ids), labels = lab.ids, 
                  xlab = "", las = las, ...)
                labline(h = 100/x$k, lty = "dashed", ...)
                if (na.act == "na.exclude" || na.act == "na.pass") 
                  llines(seq_len(len.ids)[not.na], zi[not.na], 
                    col = col.na, ...)
                llines(seq_len(len.ids), zi, ...)
                lpoints(seq_len(len.ids), zi, pch = pch, bg = bg, 
                  ...)
                lpoints(seq_len(len.ids)[ids.infl], zi[ids.infl], 
                  bg = bg.infl, pch = pch, ...)
            }
        }
    }
    if (plotdfb) {
        if (dfbnew) {
            dev.new()
            par.mar <- par("mar")
            par.mar.adj <- par.mar - c(2, 2, 2, 1)
            par.mar.adj[par.mar.adj < 1] <- 1
            par(mar = par.mar.adj)
            on.exit(par(mar = par.mar), add = TRUE)
        }
        else {
            if (plotinf) {
                par.ask <- par("ask")
                par(ask = TRUE)
                on.exit(par(ask = par.ask), add = TRUE)
            }
        }
        par(mfrow = c(length(which.dfb), 1))
        for (i in seq_len(length(which.dfb))) {
            zi <- x$dfb[, which.dfb[i]]
            not.na <- !is.na(zi)
            if (na.act == "na.omit") {
                zi <- zi[not.na]
                len.ids <- length(x$ids) - sum(!not.na)
                ids.infl <- x$is.infl[not.na]
                lab.ids <- ids[not.na]
            }
            if (na.act == "na.exclude" || na.act == "na.pass") {
                len.ids <- length(x$ids)
                ids.infl <- x$is.infl
                lab.ids <- ids
            }
            lplot(NA, NA, xlim = c(1, len.ids), ylim = range(zi, 
                na.rm = TRUE), xaxt = "n", main = paste("dfb: ", 
                colnames(x$dfb)[which.dfb[i]]), xlab = "", ylab = "", 
                las = las, ...)
            laxis(side = 1, at = seq_len(len.ids), labels = lab.ids, 
                xlab = "", las = las, ...)
            labline(h = 0, lty = "dashed", ...)
            labline(h = 1, lty = "dotted", ...)
            labline(h = -1, lty = "dotted", ...)
            if (na.act == "na.exclude" || na.act == "na.pass") 
                llines(seq_len(len.ids)[not.na], zi[not.na], 
                  col = col.na, ...)
            llines(seq_len(len.ids), zi, ...)
            lpoints(seq_len(len.ids), zi, pch = pch, bg = bg, 
                ...)
            lpoints(seq_len(len.ids)[ids.infl], zi[ids.infl], 
                bg = bg.infl, pch = pch, ...)
        }
    }
    invisible()
}
