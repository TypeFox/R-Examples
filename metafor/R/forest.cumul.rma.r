forest.cumul.rma <-
function (x, annotate = TRUE, xlim, alim, clim, ylim, at, steps = 5, 
    level = x$level, digits = 2, refline = 0, xlab, ilab, ilab.xpos, 
    ilab.pos, transf, atransf, targs, rows, efac = 1, pch = 15, 
    psize = 1, lty, cex, cex.lab, cex.axis, ...) 
{
    na.act <- getOption("na.action")
    if (!is.element("cumul.rma", class(x))) 
        stop("Argument 'x' must be an object of class \"cumul.rma\".")
    if (missing(transf)) 
        transf <- FALSE
    if (missing(atransf)) 
        atransf <- FALSE
    transf.char <- deparse(substitute(transf))
    atransf.char <- deparse(substitute(atransf))
    if (transf.char != "FALSE" && atransf.char != "FALSE") 
        stop("Use either 'transf' or 'atransf' to specify a transformation (not both).")
    if (missing(targs)) 
        targs <- NULL
    if (missing(at)) 
        at <- NULL
    if (missing(ilab)) 
        ilab <- NULL
    if (missing(ilab.xpos)) 
        ilab.xpos <- NULL
    if (missing(ilab.pos)) 
        ilab.pos <- NULL
    if (missing(cex)) 
        cex <- NULL
    if (missing(cex.lab)) 
        cex.lab <- NULL
    if (missing(cex.axis)) 
        cex.axis <- NULL
    if (missing(lty)) {
        lty <- c("solid", "solid")
    }
    else {
        if (length(lty) == 1L) 
            lty <- c(lty, "solid")
    }
    if (length(digits) == 1L) 
        digits <- c(digits, digits)
    alpha <- ifelse(level > 1, (100 - level)/100, 1 - level)
    yi <- x$estimate
    if (is.null(attr(yi, "measure"))) {
        measure <- "GEN"
    }
    else {
        measure <- attr(yi, "measure")
    }
    vi <- x$se^2
    ci.lb <- x$ci.lb
    ci.ub <- x$ci.ub
    if (length(yi) != length(vi)) 
        stop("Length of yi and vi (or sei) vectors are not the same.")
    k <- length(yi)
    if (x$slab.null) {
        slab <- paste("+ Study ", x$slab)
        slab[1] <- paste("Study ", x$slab[1])
    }
    else {
        slab <- paste("+", x$slab)
        slab[1] <- paste(x$slab[1])
    }
    if (is.vector(ilab) || NCOL(ilab) == 1) 
        ilab <- cbind(ilab)
    if (length(pch) == 1L) 
        pch <- rep(pch, k)
    if (length(pch) != length(yi)) 
        stop("Number of outcomes does not correspond to the length of the 'pch' argument.")
    if (!is.null(psize)) {
        if (length(psize) == 1L) 
            psize <- rep(psize, k)
        if (length(psize) != length(yi)) 
            stop("Number of outcomes does not correspond to the length of the 'psize' argument.")
    }
    if (missing(rows)) {
        rows <- k:1
    }
    else {
        if (length(rows) == 1L) 
            rows <- rows:(rows - k + 1)
    }
    if (length(rows) != length(yi)) 
        stop("Number of outcomes does not correspond to the length of the 'rows' argument.")
    yi <- yi[k:1]
    vi <- vi[k:1]
    ci.lb <- ci.lb[k:1]
    ci.ub <- ci.ub[k:1]
    slab <- slab[k:1]
    ilab <- ilab[k:1, , drop = FALSE]
    pch <- pch[k:1]
    psize <- psize[k:1]
    rows <- rows[k:1]
    yivi.na <- is.na(yi) | is.na(vi)
    if (any(yivi.na)) {
        not.na <- !yivi.na
        if (na.act == "na.omit") {
            yi <- yi[not.na]
            vi <- vi[not.na]
            ci.lb <- ci.lb[not.na]
            ci.ub <- ci.ub[not.na]
            slab <- slab[not.na]
            ilab <- ilab[not.na, , drop = FALSE]
            pch <- pch[not.na]
            psize <- psize[not.na]
            rows.new <- rows
            rows.na <- rows[!not.na]
            for (j in seq_len(length(rows.na))) {
                rows.new[rows >= rows.na[j]] <- rows.new[rows >= 
                  rows.na[j]] - 1
            }
            rows <- rows.new[not.na]
        }
        if (na.act == "na.fail") 
            stop("Missing values in results.")
    }
    k <- length(yi)
    if (is.function(transf)) {
        if (is.null(targs)) {
            yi <- sapply(yi, transf)
            ci.lb <- sapply(ci.lb, transf)
            ci.ub <- sapply(ci.ub, transf)
        }
        else {
            yi <- sapply(yi, transf, targs)
            ci.lb <- sapply(ci.lb, transf, targs)
            ci.ub <- sapply(ci.ub, transf, targs)
        }
    }
    ci.bounds <- cbind(ci.lb, ci.ub)
    rev.order <- ifelse(ci.ub < ci.lb, TRUE, FALSE)
    rev.order[is.na(rev.order)] <- FALSE
    ci.bounds[rev.order, ] <- ci.bounds[rev.order, 2:1]
    ci.lb <- ci.bounds[, 1]
    ci.ub <- ci.bounds[, 2]
    if (!missing(clim)) {
        clim <- sort(clim)
        if (length(clim) != 2L) 
            stop("Argument 'clim' must be of length 2.")
        ci.lb[ci.lb < clim[1]] <- clim[1]
        ci.ub[ci.ub > clim[2]] <- clim[2]
    }
    if (is.null(psize)) {
        if (any(vi <= 0, na.rm = TRUE)) {
            psize <- rep(1, k)
        }
        else {
            wi <- 1/sqrt(vi)
            psize <- wi/sum(wi, na.rm = TRUE)
            psize <- (psize - min(psize, na.rm = TRUE))/(max(psize, 
                na.rm = TRUE) - min(psize, na.rm = TRUE))
            psize <- (psize * 1) + 0.5
            if (all(is.na(psize))) 
                psize <- rep(1, k)
        }
    }
    rng <- max(ci.ub, na.rm = TRUE) - min(ci.lb, na.rm = TRUE)
    if (annotate) {
        plot.multp.l <- 1.2
        plot.multp.r <- 1.2
        axis.multp.l <- 0.2
        axis.multp.r <- 0.2
    }
    else {
        plot.multp.l <- 1.2
        plot.multp.r <- 0.4
        axis.multp.l <- 0.2
        axis.multp.r <- 0.2
    }
    if (missing(xlim)) {
        xlim <- c(min(ci.lb, na.rm = TRUE) - rng * plot.multp.l, 
            max(ci.ub, na.rm = TRUE) + rng * plot.multp.r)
        xlim <- round(xlim, digits[2])
    }
    alim.spec <- TRUE
    if (missing(alim)) {
        if (is.null(at)) {
            alim <- range(pretty(x = c(min(ci.lb, na.rm = TRUE), 
                max(ci.ub, na.rm = TRUE)), n = steps - 1))
            alim.spec <- FALSE
        }
        else {
            alim <- range(at)
        }
    }
    alim <- sort(alim)
    xlim <- sort(xlim)
    if (xlim[1] > min(yi, na.rm = TRUE)) {
        xlim[1] <- min(yi, na.rm = TRUE)
    }
    if (xlim[2] < max(yi, na.rm = TRUE)) {
        xlim[2] <- max(yi, na.rm = TRUE)
    }
    if (alim[1] < xlim[1]) {
        xlim[1] <- alim[1]
    }
    if (alim[2] > xlim[2]) {
        xlim[2] <- alim[2]
    }
    if (missing(ylim)) {
        ylim <- c(0.5, k + 3)
    }
    else {
        ylim <- sort(ylim)
    }
    if (is.null(at)) {
        if (alim.spec) {
            at <- seq(from = alim[1], to = alim[2], length.out = steps)
        }
        else {
            at <- pretty(x = c(min(ci.lb, na.rm = TRUE), max(ci.ub, 
                na.rm = TRUE)), n = steps - 1)
        }
    }
    else {
        at[at < alim[1]] <- alim[1]
        at[at > alim[2]] <- alim[2]
        at <- unique(at)
    }
    at.lab <- at
    if (is.function(atransf)) {
        if (is.null(targs)) {
            at.lab <- formatC(sapply(at.lab, atransf), digits = digits[2], 
                format = "f")
        }
        else {
            at.lab <- formatC(sapply(at.lab, atransf, targs), 
                digits = digits[2], format = "f")
        }
    }
    else {
        at.lab <- formatC(at.lab, digits = digits[2], format = "f")
    }
    par.mar <- par("mar")
    par.mar.adj <- par.mar - c(0, 3, 1, 1)
    par.mar.adj[par.mar.adj < 0] <- 0
    par(mar = par.mar.adj)
    on.exit(par(mar = par.mar))
    plot(NA, NA, xlim = xlim, ylim = ylim, xlab = "", ylab = "", 
        yaxt = "n", xaxt = "n", xaxs = "i", bty = "n", ...)
    abline(h = ylim[2] - 2, lty = lty[2], ...)
    if (is.numeric(refline)) 
        segments(refline, ylim[1] - 5, refline, ylim[2] - 2, 
            lty = "dotted", ...)
    par.usr <- par("usr")
    height <- par.usr[4] - par.usr[3]
    lheight <- strheight("O")
    cex.adj <- ifelse(k * lheight > height * 0.8, height/(1.25 * 
        k * lheight), 1)
    if (is.null(cex)) {
        cex <- par("cex") * cex.adj
    }
    else {
        if (is.null(cex.lab)) 
            cex.lab <- cex
        if (is.null(cex.axis)) 
            cex.axis <- cex
    }
    if (is.null(cex.lab)) 
        cex.lab <- par("cex") * cex.adj
    if (is.null(cex.axis)) 
        cex.axis <- par("cex") * cex.adj
    axis(side = 1, at = at, labels = at.lab, cex.axis = cex.axis, 
        ...)
    if (missing(xlab)) 
        xlab <- .setlab(measure, transf.char, atransf.char, gentype = 2)
    mtext(xlab, side = 1, at = min(at) + (max(at) - min(at))/2, 
        line = par("mgp")[1] - 0.5, cex = cex.lab, ...)
    for (i in seq_len(k)) {
        if (is.na(yi[i]) || is.na(vi[i])) 
            next
        if (ci.lb[i] >= alim[2]) {
            polygon(x = c(alim[2], alim[2] - (1.4/100) * cex * 
                (xlim[2] - xlim[1]), alim[2] - (1.4/100) * cex * 
                (xlim[2] - xlim[1]), alim[2]), y = c(rows[i], 
                rows[i] + (height/150) * cex * efac, rows[i] - 
                  (height/150) * cex * efac, rows[i]), col = "black", 
                ...)
            next
        }
        if (ci.ub[i] <= alim[1]) {
            polygon(x = c(alim[1], alim[1] + (1.4/100) * cex * 
                (xlim[2] - xlim[1]), alim[1] + (1.4/100) * cex * 
                (xlim[2] - xlim[1]), alim[1]), y = c(rows[i], 
                rows[i] + (height/150) * cex * efac, rows[i] - 
                  (height/150) * cex * efac, rows[i]), col = "black", 
                ...)
            next
        }
        segments(max(ci.lb[i], alim[1]), rows[i], min(ci.ub[i], 
            alim[2]), rows[i], lty = lty[1], ...)
        if (ci.lb[i] >= alim[1]) {
            segments(ci.lb[i], rows[i] - (height/150) * cex * 
                efac, ci.lb[i], rows[i] + (height/150) * cex * 
                efac, ...)
        }
        else {
            polygon(x = c(alim[1], alim[1] + (1.4/100) * cex * 
                (xlim[2] - xlim[1]), alim[1] + (1.4/100) * cex * 
                (xlim[2] - xlim[1]), alim[1]), y = c(rows[i], 
                rows[i] + (height/150) * cex * efac, rows[i] - 
                  (height/150) * cex * efac, rows[i]), col = "black", 
                ...)
        }
        if (ci.ub[i] <= alim[2]) {
            segments(ci.ub[i], rows[i] - (height/150) * cex * 
                efac, ci.ub[i], rows[i] + (height/150) * cex * 
                efac, ...)
        }
        else {
            polygon(x = c(alim[2], alim[2] - (1.4/100) * cex * 
                (xlim[2] - xlim[1]), alim[2] - (1.4/100) * cex * 
                (xlim[2] - xlim[1]), alim[2]), y = c(rows[i], 
                rows[i] + (height/150) * cex * efac, rows[i] - 
                  (height/150) * cex * efac, rows[i]), col = "black", 
                ...)
        }
    }
    text(xlim[1], rows, slab, pos = 4, cex = cex, ...)
    if (!is.null(ilab)) {
        if (is.null(ilab.xpos)) 
            stop("Must specify 'ilab.xpos' argument when adding information with 'ilab'.")
        if (length(ilab.xpos) != NCOL(ilab)) 
            stop("Number of 'ilab' columns does not match length of 'ilab.xpos' argument.")
        for (l in seq_len(NCOL(ilab))) {
            text(ilab.xpos[l], rows, ilab[, l], pos = ilab.pos[l], 
                cex = cex, ...)
        }
    }
    if (annotate) {
        if (is.function(atransf)) {
            if (is.null(targs)) {
                annotext <- round(cbind(sapply(yi, atransf), 
                  sapply(ci.lb, atransf), sapply(ci.ub, atransf)), 
                  digits[1])
            }
            else {
                annotext <- round(cbind(sapply(yi, atransf, targs), 
                  sapply(ci.lb, atransf, targs), sapply(ci.ub, 
                    atransf, targs)), digits[1])
            }
            rev.order <- ifelse(annotext[, 3] < annotext[, 2], 
                TRUE, FALSE)
            rev.order[is.na(rev.order)] <- FALSE
            annotext[rev.order, 2:3] <- annotext[rev.order, 3:2]
        }
        else {
            annotext <- round(cbind(yi, ci.lb, ci.ub), digits[1])
        }
        annotext <- matrix(apply(annotext, 2, format, nsmall = digits[1]), 
            ncol = 3)
        annotext <- cbind(annotext[, 1], " [ ", annotext[, 2], 
            " , ", annotext[, 3], " ]")
        annotext <- apply(annotext, 1, paste, collapse = "")
        text(x = xlim[2], rows, labels = annotext, pos = 2, cex = cex, 
            ...)
    }
    for (i in seq_len(k)) {
        if (is.na(yi[i])) 
            next
        if (yi[i] >= alim[1] && yi[i] <= alim[2]) 
            points(yi[i], rows[i], pch = pch[i], cex = cex * 
                psize[i], ...)
    }
    res <- list(xlim = par("usr")[1:2], alim = alim, at = at, 
        ylim = ylim, rows = rows, cex = cex, cex.lab = cex.lab, 
        cex.axis = cex.axis)
    invisible(res)
}
