plot.orlm <- function (x, caption = "Residuals vs Fitted", 
    panel = if (add.smooth) panel.smooth else points, 
    sub.caption = NULL, main = "", ..., id.n = 3, labels.id = names(x$residuals), 
    cex.id = 0.75, add.smooth = getOption("add.smooth"), 
    label.pos = c(4, 2), cex.caption = 1) 
{
    ## adapted from plot.lm
    if (!inherits(x, "orlm")) 
        stop("use only with \"orlm\" objects")
    if (is.null(x$residuals)) 
        stop("plot only works on orlm objects that contain residuals and fitted values")
    isGlm <- inherits(x, "orglm")
    r <- x$residuals
    yh <- x$fitted.values
    w <- x$weights
    if (!is.null(w)) {
        wind <- w != 0
        r <- r[wind]
        yh <- yh[wind]
        w <- w[wind]
        labels.id <- labels.id[wind]
    }
    n <- length(r)
        l.fit <- if (isGlm) 
            "Predicted values"
        else "Fitted values"
    if (is.null(id.n)) 
        id.n <- 0
    else {
        id.n <- as.integer(id.n)
        if (id.n < 0 || id.n > n) 
            stop(gettextf("'id.n' must be in {1,..,%d}", n), 
                domain = NA)
    }
    if (id.n > 0) {
        if (is.null(labels.id)) 
            labels.id <- paste(1:n)
        iid <- 1:id.n
        show.r <- sort.list(abs(r), decreasing = TRUE)[iid]
        text.id <- function(x, y, ind, adj.x = TRUE) {
            labpos <- if (adj.x) 
                label.pos[1 + as.numeric(x > mean(range(x)))]
            else 3
            text(x, y, labels.id[ind], cex = cex.id, xpd = TRUE, 
                pos = labpos, offset = 0.25)
        }
    }
    if (is.null(sub.caption)) {
        cal <- "Order-restricted linear model"
        if (x$meq>0) cal <- paste(cal, "with", x$meq, "equality and", nrow(x$ui)-x$meq, 
            "inequality restrictions of which", length(x$iact)-x$meq, "are active.") 
        else cal <- paste(cal, "with", nrow(x$ui)-x$meq, 
            "inequality restrictions of which", length(x$iact), "are active.")
        sub.caption <- cal 
    }
    one.fig <- prod(par("mfcol")) == 1
        ylim <- range(r, na.rm = TRUE)
        if (id.n > 0) 
            ylim <- extendrange(r = ylim, f = 0.08)
        plot(yh, r, xlab = l.fit, ylab = "Residuals", main = main, 
            ylim = ylim, type = "n", ...)
        panel(yh, r, ...)
        if (one.fig) 
            title(sub = sub.caption, ...)
        mtext(caption[1], 3, 0.25, cex = cex.caption)
        if (id.n > 0) {
            y.id <- r[show.r]
            y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
            text.id(yh[show.r], y.id, show.r)
        }
        abline(h = 0, lty = 3, col = "gray")

   if (!one.fig && par("oma")[3] >= 1) 
        mtext(sub.caption, outer = TRUE, cex = 1.25)
    invisible()
}
