plot.hdlm <- function (x, which = c(1L:4L), caption = list("Residuals vs Fitted", 
    "Normal Q-Q", "Distribution of Coefficients", "P-value vs Confidence Intervals"), 
    panel = if (add.smooth) panel.smooth else points, sub.caption = NULL, 
    main = "", ask = prod(par("mfcol")) < length(which) && dev.interactive(), 
    ..., id.n = 3, labels.id = names(residuals(x)), cex.id = 0.75, 
    qqline = TRUE, cook.levels = c(0.5, 1), add.smooth = getOption("add.smooth"), 
    label.pos = c(4, 2), cex.caption = 1) 
{
    dropInf <- function(x, h) {
        if (any(isInf <- h >= 1)) {
            warning("Not plotting observations with leverage one:\n  ", 
                paste(which(isInf), collapse = ", "), call. = FALSE)
            x[isInf] <- NaN
        }
        x
    }
    if (!inherits(x, "hdlm")) 
        stop("use only with \"hdlm\" objects")
    if (!is.numeric(which) || any(which < 1) || any(which > 6)) 
        stop("'which' must be in 1:6")
    show <- rep(FALSE, 6)
    show[which] <- TRUE
    r <- residuals(x)
    yh <- predict(x)
    n <- length(r)
    l.fit <- "Fitted values"
    if (is.null(id.n)) 
        id.n <- 0
    else {
        id.n <- as.integer(id.n)
        if (id.n < 0L || id.n > n) 
            stop(gettextf("'id.n' must be in {1,..,%d}", n), 
                domain = NA)
    }
    if (id.n > 0L) {
        if (is.null(labels.id)) 
            labels.id <- paste(1L:n)
        iid <- 1L:id.n
        show.r <- sort.list(abs(r), decreasing = TRUE)[iid]
        text.id <- function(x, y, ind, adj.x = TRUE) {
            labpos <- if (adj.x) 
                label.pos[1 + as.numeric(x > mean(range(x)))]
            else 3
            text(x, y, labels.id[ind], cex = cex.id, xpd = TRUE, 
                pos = labpos, offset = 0.25)
        }
    }
    getCaption <- function(k) if (length(caption) < k) 
        NA_character_
    else as.graphicsAnnot(caption[[k]])
    if (is.null(sub.caption)) {
        cal <- x$call
        if (!is.na(m.f <- match("formula", names(cal)))) {
            cal <- cal[c(1, m.f)]
            names(cal)[2L] <- ""
        }
        cc <- deparse(cal, 80)
        nc <- nchar(cc[1L], "c")
        abbr <- length(cc) > 1 || nc > 75
        sub.caption <- if (abbr) 
            paste(substr(cc[1L], 1L, min(75L, nc)), "...")
        else cc[1L]
    }
    one.fig <- prod(par("mfcol")) == 1
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    if (show[1L]) {
        ylim <- range(r, na.rm = TRUE)
        if (id.n > 0) 
            ylim <- extendrange(r = ylim, f = 0.08)
        dev.hold()
        plot(yh, r, xlab = l.fit, ylab = "Residuals", main = main, 
            ylim = ylim, type = "n", ...)
        panel(yh, r, ...)
        if (one.fig) 
            title(sub = sub.caption, ...)
        mtext(getCaption(1), 3, 0.25, cex = cex.caption)
        if (id.n > 0) {
            y.id <- r[show.r]
            y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
            text.id(yh[show.r], y.id, show.r)
        }
        abline(h = 0, lty = 3, col = "gray")
        dev.flush()
    }
    if (show[2L]) {
        ylim <- range(r, na.rm = TRUE)
        ylim[2L] <- ylim[2L] + diff(ylim) * 0.075
        dev.hold()
        qq <- qqnorm(r, main = main, ylab = "Residuals", ylim = ylim, 
            ...)
        if (qqline) 
            qqline(r, lty = 3, col = "gray50")
        if (one.fig) 
            title(sub = sub.caption, ...)
        mtext(getCaption(2), 3, 0.25, cex = cex.caption)
        if (id.n > 0) 
            text.id(qq$x[show.r], qq$y[show.r], show.r)
        dev.flush()
    }
    if (show[3L]) {
        DIFF <- (x$upper.bound - x$lower.bound)
        vals <- x$coef / DIFF
        vals[!is.finite(vals)] <- 0
        dev.hold()
        ylim <- range(vals, na.rm = TRUE)
        ylim[1L] <- ylim[1L] - diff(ylim) * 0.075
        ylim[2L] <- ylim[2L] + diff(ylim) * 0.075
        qq <- barplot(sort(vals), col="gray75", names.arg="", ylab="Standardized Coefficients",
                        main = main, ylim=ylim, xlab="Sorted Coefficients", ...)
        mtext(getCaption(3), 3, 0.25, cex = cex.caption)
        box()
        dev.flush()
    }
    if (show[4L]) {
        DIFF <- (x$upper.bound - x$lower.bound)
        DIFF <- DIFF / 2
        SE <- DIFF / qnorm(1-x$siglevel/2)
        pval_new <- (1 - pnorm( abs(x$coef / SE) ))*2
        DIFF <- (x$upper.bound - x$lower.bound)
        ylim <- range(pval_new, na.rm = TRUE)
        ylim[2L] <- ylim[2L] + diff(ylim) * 0.075
        pval_new[!is.finite(pval_new)] <- 1
        dev.hold()
        plot(x$p.value, pval_new, ylab = "Normal Model P-value", xlab="Reported P-value",
            main = main, ylim = ylim, type = "p",  ...)
        if (one.fig) 
            title(sub = sub.caption, ...)
        mtext(getCaption(4), 3, 0.25, cex = cex.caption)
        abline(0,1,col='red')
        if (id.n > 0) {
            y.id <- pval_new[show.r]
            y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
            text.id(x$p.value[show.r], y.id, show.r)
        }
        dev.flush()
    }
    if (!one.fig && par("oma")[3L] >= 1) 
        mtext(sub.caption, outer = TRUE, cex = 1.25)
    invisible()
}

