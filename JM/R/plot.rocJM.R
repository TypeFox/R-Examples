plot.rocJM <-
function (x, which = NULL, type = c("ROC", "AUC"), 
    ndt = "all", main = NULL, caption = NULL, xlab = NULL, ylab = NULL, 
    ask = NULL, legend = FALSE, lx = NULL, ly = NULL, lty = NULL, col = NULL, 
    cex.caption = 0.8, cex.axis = NULL, cex.lab = NULL, cex.main = NULL, ...) {
    type <- match.arg(type)
    if (is.null(which))
        which <- grep("MCresults", names(x))
    if (is.null(ask))
        ask <- prod(par("mfcol")) < length(which)
    if (ask) {
        op <- par(ask = TRUE)
        on.exit(par(op))
    }
    nams <- if (is.matrix(aucs <- x$AUCs)) colnames(aucs) else names(aucs)
    rnams <- if (is.matrix(aucs)) rownames(aucs) else lapply(aucs, names)
    if (is.null(main))
        main <- paste("Case", names(x$times))
    if (is.null(caption))
        caption <- sapply(x$times, function (xx) {
            paste("(Follow-up time(s): ", 
                paste(round(xx, 1), collapse = ", "), ")", sep = "")
        })
    if (type == "ROC") {
        if (is.null(xlab))
            xlab <- "1 - Specificity"
        if (is.null(ylab))
            ylab <- "Sensitivity"
        if (is.null(lx))
            lx <- "bottomright"
        if (is.null(lty))
            lty <- 1
        for (i in seq_along(which)) {
            ii <- which[i]
            rr <- x[[ii]]
            col. <- if (is.null(col)) rep(1:6, length.out = nrow(rr$Sens)) else col
            if (ndt == "all") {
                matplot(t(1 - rr$Spec), t(rr$Sen), type = "l", lty = lty, col = col., 
                    ylim = c(0, 1), xlim = c(0, 1), xlab = xlab, 
                    ylab = ylab, main = main[ii], cex.axis = cex.axis, 
                    cex.lab = cex.lab, cex.main = cex.main, ...)
                abline(a = 0, b = 1, col = "grey")
                mtext(caption[ii], 3, 0.25, cex = cex.caption)
                if (legend) {
                    labs <- if (is.list(rnams)) rnams[[i]] else rnams
                    labs <- sapply(rnams, function (n)
                        as.expression(substitute(paste(Delta, "t = ", n), list("n" = n))))
                    legend(lx, ly, labs, lty = lty, col = col., bty = "n", ...)
                }
            } else {
                ind <- if (ndt > (nr <- nrow(rr$Sen))) 1:nr else 
                    round(seq(1, nr, length.out = ndt))
                matplot(t(1 - rr$Spec)[, ind], t(rr$Sen)[, ind], 
                    type = "l", lty = lty, col = col., ylim = c(0, 1), 
                    xlim = c(0, 1), xlab = xlab, ylab = ylab, main = main[ii], 
                    cex.axis = cex.axis, cex.lab = cex.lab, cex.main = cex.main, ...)
                abline(a = 0, b = 1, col = "grey")            
                mtext(caption[ii], 3, 0.25, cex = cex.caption)
                if (legend) {
                    labs <- if (is.list(rnams)) rnams[[i]] else rnams 
                    labs <- sapply(rnams, function (n)
                        as.expression(substitute(paste(Delta, "t = ", n), list("n" = n))))
                    legend(lx, ly, labs[ind], lty = lty, col = col., bty = "n", ...)
                }
            } 
        }
    } else {
        if (is.null(xlab))
            xlab <- "Time"
        if (is.null(ylab))
            ylab <- "AUC"
        if (is.null(lx))
            lx <- "topleft"
        if (is.null(lty))
            lty <- 1
        if (is.list(aucs)) {
            max.n <- max(sapply(aucs, length))
            aucs <- sapply(aucs, function (y) {
                out <- rep(as.numeric(NA), max.n)
                names(out) <- sprintf("%.1f", 
                    unique(sort(as.numeric(unlist(rnams)))))
                out[match(names(y), names(out))] <- y
                out
            })
        }
        col. <- if (is.null(col)) rep(1:6, length.out = ncol(aucs)) else col
        ind.which <- if (is.numeric(which)) which else match(which, nams)
        nr <- nrow(aucs)
        ind.time <- if (ndt == "all" || ndt > nr) 1:nr else 
            round(seq(1, nr, length.out = ndt))
        matplot(as.numeric(rownames(aucs))[ind.time], 
            aucs[ind.time, ind.which, drop = FALSE], type = "l", 
            lty = lty, col = col., xlab = xlab, ylab = ylab, 
            cex.axis = cex.axis, cex.lab = cex.lab, ...)
        if (legend)
            legend(lx, ly, main, lty = lty, col = col., bty = "n", ...)
    }
    invisible()
}
