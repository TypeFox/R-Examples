forest.rma <-
function (x, annotate = TRUE, addfit = TRUE, addcred = FALSE, 
    showweights = FALSE, xlim, alim, clim, ylim, at, steps = 5, 
    level = x$level, digits = 2, refline = 0, xlab, slab, mlab, 
    ilab, ilab.xpos, ilab.pos, order, transf, atransf, targs, 
    rows, efac = 1, pch = 15, psize, col, border, lty, cex, cex.lab, 
    cex.axis, ...) 
{
    if (!is.element("rma", class(x))) 
        stop("Argument 'x' must be an object of class \"rma\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
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
    if (missing(order)) 
        order <- NULL
    if (missing(psize)) 
        psize <- NULL
    if (missing(cex)) 
        cex <- NULL
    if (missing(cex.lab)) 
        cex.lab <- NULL
    if (missing(cex.axis)) 
        cex.axis <- NULL
    if (x$int.only) {
        if (missing(col)) {
            col <- c("black", "gray50")
        }
        else {
            if (length(col) == 1L) 
                col <- c(col, "gray50")
        }
        if (missing(border)) 
            border <- "black"
    }
    else {
        if (missing(col)) 
            col <- "gray"
        if (missing(border)) 
            border <- "gray"
    }
    if (missing(lty)) {
        lty <- c("solid", "dotted", "solid")
    }
    else {
        if (length(lty) == 1L) 
            lty <- c(lty, "dotted", "solid")
        if (length(lty) == 2L) 
            lty <- c(lty, "solid")
    }
    if (length(efac) == 1) 
        efac <- rep(efac, 2)
    measure <- x$measure
    if (is.element("rma.glmm", class(x)) && showweights) 
        stop("Option 'showweights=TRUE' currently not possible for 'rma.glmm' objects. Sorry!")
    if (length(digits) == 1L) 
        digits <- c(digits, digits)
    alpha <- ifelse(level > 1, (100 - level)/100, 1 - level)
    yi <- x$yi.f
    vi <- x$vi.f
    X <- x$X.f
    k <- length(yi)
    if (missing(slab)) {
        if (x$slab.null) {
            slab <- paste("Study ", x$slab)
        }
        else {
            slab <- x$slab
        }
    }
    if (length(yi) != length(slab)) 
        stop("Number of outcomes does not correspond to the length of the 'slab' argument.")
    if (is.vector(ilab) || NCOL(ilab) == 1) 
        ilab <- cbind(ilab)
    if (length(pch) == 1L) 
        pch <- rep(pch, k)
    if (length(pch) != length(yi)) 
        stop("Number of outcomes does not correspond to the length of the 'pch' argument.")
    options(na.action = "na.pass")
    if (x$int.only) {
        pred <- fitted(x)
        pred.ci.lb <- rep(NA_real_, k)
        pred.ci.ub <- rep(NA_real_, k)
    }
    else {
        temp <- predict(x, level = level)
        pred <- temp$pred
        if (addcred) {
            pred.ci.lb <- temp$cr.lb
            pred.ci.ub <- temp$cr.ub
        }
        else {
            pred.ci.lb <- temp$ci.lb
            pred.ci.ub <- temp$ci.ub
        }
    }
    if (is.element("rma.glmm", class(x))) {
        weights <- NULL
    }
    else {
        weights <- weights(x)
    }
    options(na.action = na.act)
    if (!is.null(psize)) {
        if (length(psize) == 1L) 
            psize <- rep(psize, k)
        if (length(psize) != length(yi)) 
            stop("Number of outcomes does not correspond to the length of the 'psize' argument.")
    }
    if (!is.null(order)) {
        if (is.character(order)) {
            if (length(order) != 1) 
                stop("Incorrect length of order argument.")
            if (order == "obs") 
                sort.vec <- order(yi)
            if (order == "fit") 
                sort.vec <- order(pred)
            if (order == "prec") 
                sort.vec <- order(vi, yi)
            if (order == "resid") 
                sort.vec <- order(yi - pred, yi)
            if (order == "rstandard") 
                sort.vec <- order(rstandard(x)$z, yi)
            if (order == "abs.resid") 
                sort.vec <- order(abs(yi - pred), yi)
            if (order == "abs.rstandard") 
                sort.vec <- order(abs(rstandard(x)$z), yi)
        }
        else {
            sort.vec <- order
        }
        yi <- yi[sort.vec]
        vi <- vi[sort.vec]
        X <- X[sort.vec, , drop = FALSE]
        slab <- slab[sort.vec]
        ilab <- ilab[sort.vec, , drop = FALSE]
        pred <- pred[sort.vec]
        pred.ci.lb <- pred.ci.lb[sort.vec]
        pred.ci.ub <- pred.ci.ub[sort.vec]
        weights <- weights[sort.vec]
        pch <- pch[sort.vec]
        psize <- psize[sort.vec]
    }
    k <- length(yi)
    if (missing(rows)) {
        rows <- k:1
    }
    else {
        if (length(rows) == 1L) {
            rows <- rows:(rows - k + 1)
        }
    }
    if (length(rows) != length(yi)) 
        stop("Number of outcomes does not correspond to the length of the 'rows' argument.")
    yi <- yi[k:1]
    vi <- vi[k:1]
    X <- X[k:1, , drop = FALSE]
    slab <- slab[k:1]
    ilab <- ilab[k:1, , drop = FALSE]
    pred <- pred[k:1]
    pred.ci.lb <- pred.ci.lb[k:1]
    pred.ci.ub <- pred.ci.ub[k:1]
    weights <- weights[k:1]
    pch <- pch[k:1]
    psize <- psize[k:1]
    rows <- rows[k:1]
    yiviX.na <- is.na(yi) | is.na(vi) | apply(is.na(X), 1, any)
    if (any(yiviX.na)) {
        not.na <- !yiviX.na
        if (na.act == "na.omit") {
            yi <- yi[not.na]
            vi <- vi[not.na]
            X <- X[not.na, , drop = FALSE]
            slab <- slab[not.na]
            ilab <- ilab[not.na, , drop = FALSE]
            pred <- pred[not.na]
            pred.ci.lb <- pred.ci.lb[not.na]
            pred.ci.ub <- pred.ci.ub[not.na]
            weights <- weights[not.na]
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
    ci.lb <- yi - qnorm(alpha/2, lower.tail = FALSE) * sqrt(vi)
    ci.ub <- yi + qnorm(alpha/2, lower.tail = FALSE) * sqrt(vi)
    if (is.function(transf)) {
        if (is.null(targs)) {
            yi <- sapply(yi, transf)
            ci.lb <- sapply(ci.lb, transf)
            ci.ub <- sapply(ci.ub, transf)
            pred <- sapply(pred, transf)
            pred.ci.lb <- sapply(pred.ci.lb, transf)
            pred.ci.ub <- sapply(pred.ci.ub, transf)
        }
        else {
            yi <- sapply(yi, transf, targs)
            ci.lb <- sapply(ci.lb, transf, targs)
            ci.ub <- sapply(ci.ub, transf, targs)
            pred <- sapply(pred, transf, targs)
            pred.ci.lb <- sapply(pred.ci.lb, transf, targs)
            pred.ci.ub <- sapply(pred.ci.ub, transf, targs)
        }
    }
    ci.bounds <- cbind(ci.lb, ci.ub)
    rev.order <- ifelse(ci.ub < ci.lb, TRUE, FALSE)
    rev.order[is.na(rev.order)] <- FALSE
    ci.bounds[rev.order, ] <- ci.bounds[rev.order, 2:1]
    ci.lb <- ci.bounds[, 1]
    ci.ub <- ci.bounds[, 2]
    pred.ci.bounds <- cbind(pred.ci.lb, pred.ci.ub)
    rev.order <- ifelse(pred.ci.ub < pred.ci.lb, TRUE, FALSE)
    rev.order[is.na(rev.order)] <- FALSE
    pred.ci.bounds[rev.order, ] <- pred.ci.bounds[rev.order, 
        2:1]
    pred.ci.lb <- pred.ci.bounds[, 1]
    pred.ci.ub <- pred.ci.bounds[, 2]
    if (!missing(clim)) {
        clim <- sort(clim)
        if (length(clim) != 2L) 
            stop("Argument 'clim' must be of length 2.")
        ci.lb[ci.lb < clim[1]] <- clim[1]
        ci.ub[ci.ub > clim[2]] <- clim[2]
        pred.ci.lb[pred.ci.lb < clim[1]] <- clim[1]
        pred.ci.ub[pred.ci.ub > clim[2]] <- clim[2]
    }
    if (is.null(psize)) {
        if (is.null(weights)) {
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
        else {
            wi <- weights
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
        if (showweights) {
            plot.multp.l <- 2
            plot.multp.r <- 2
            axis.multp.l <- 0.2
            axis.multp.r <- 0.2
        }
        else {
            plot.multp.l <- 1.2
            plot.multp.r <- 1.2
            axis.multp.l <- 0.2
            axis.multp.r <- 0.2
        }
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
        if (x$int.only && addfit) {
            ylim <- c(-1.5, k + 3)
        }
        else {
            ylim <- c(0.5, k + 3)
        }
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
    abline(h = ylim[2] - 2, lty = lty[3], ...)
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
    if (addfit && !x$int.only) {
        for (i in seq_len(k)) {
            if (is.na(pred[i])) 
                next
            polygon(x = c(max(pred.ci.lb[i], alim[1]), pred[i], 
                min(pred.ci.ub[i], alim[2]), pred[i]), y = c(rows[i], 
                rows[i] + (height/100) * cex * efac[2], rows[i], 
                rows[i] - (height/100) * cex * efac[2]), col = col, 
                border = border, ...)
        }
    }
    if (addfit && x$int.only) {
        if (is.element("rma.mv", class(x)) && x$withG && x$tau2s > 
            1) {
            if (!is.logical(addcred)) {
                if (length(addcred) == 1) 
                  addcred <- c(addcred, addcred)
                temp <- predict(x, level = level, tau2.levels = addcred[1], 
                  gamma2.levels = addcred[2])
                addcred <- TRUE
            }
            else {
                if (addcred) {
                  stop("Need to specify the level of the inner factor(s) via the 'addcred' argument.")
                }
                else {
                  temp <- predict(x, level = level, tau2.levels = 1, 
                    gamma2.levels = 1)
                }
            }
        }
        else {
            temp <- predict(x, level = level)
        }
        b <- temp$pred
        b.ci.lb <- temp$ci.lb
        b.ci.ub <- temp$ci.ub
        b.cr.lb <- temp$cr.lb
        b.cr.ub <- temp$cr.ub
        if (is.function(transf)) {
            if (is.null(targs)) {
                b <- sapply(b, transf)
                b.ci.lb <- sapply(b.ci.lb, transf)
                b.ci.ub <- sapply(b.ci.ub, transf)
                b.cr.lb <- sapply(b.cr.lb, transf)
                b.cr.ub <- sapply(b.cr.ub, transf)
            }
            else {
                b <- sapply(b, transf, targs)
                b.ci.lb <- sapply(b.ci.lb, transf, targs)
                b.ci.ub <- sapply(b.ci.ub, transf, targs)
                b.cr.lb <- sapply(b.cr.lb, transf, targs)
                b.cr.ub <- sapply(b.cr.ub, transf, targs)
            }
        }
        b.ci.bounds <- cbind(b.ci.lb, b.ci.ub)
        rev.order <- ifelse(b.ci.ub < b.ci.lb, TRUE, FALSE)
        rev.order[is.na(rev.order)] <- FALSE
        b.ci.bounds[rev.order, ] <- b.ci.bounds[rev.order, 2:1]
        b.ci.lb <- b.ci.bounds[, 1]
        b.ci.ub <- b.ci.bounds[, 2]
        b.cr.bounds <- cbind(b.cr.lb, b.cr.ub)
        rev.order <- ifelse(b.cr.ub < b.cr.lb, TRUE, FALSE)
        rev.order[is.na(rev.order)] <- FALSE
        b.cr.bounds[rev.order, ] <- b.cr.bounds[rev.order, 2:1]
        b.cr.lb <- b.cr.bounds[, 1]
        b.cr.ub <- b.cr.bounds[, 2]
        if (!missing(clim)) {
            b.ci.lb[b.ci.lb < clim[1]] <- clim[1]
            b.ci.ub[b.ci.ub > clim[2]] <- clim[2]
            b.cr.lb[b.cr.lb < clim[1]] <- clim[1]
            b.cr.ub[b.cr.ub > clim[2]] <- clim[2]
        }
        if (x$method != "FE" && addcred) {
            segments(max(b.cr.lb, alim[1]), -1, min(b.cr.ub, 
                alim[2]), -1, lty = lty[2], col = col[2], ...)
            if (b.cr.lb >= alim[1]) {
                segments(b.cr.lb, -1 - (height/150) * cex * efac[1], 
                  b.cr.lb, -1 + (height/150) * cex * efac[1], 
                  col = col[2], ...)
            }
            else {
                polygon(x = c(alim[1], alim[1] + (1.4/100) * 
                  cex * (xlim[2] - xlim[1]), alim[1] + (1.4/100) * 
                  cex * (xlim[2] - xlim[1]), alim[1]), y = c(-1, 
                  -1 + (height/150) * cex * efac[1], -1 - (height/150) * 
                    cex * efac[1], -1), col = col[2], border = col[2], 
                  ...)
            }
            if (b.cr.ub <= alim[2]) {
                segments(b.cr.ub, -1 - (height/150) * cex * efac[1], 
                  b.cr.ub, -1 + (height/150) * cex * efac[1], 
                  col = col[2], ...)
            }
            else {
                polygon(x = c(alim[2], alim[2] - (1.4/100) * 
                  cex * (xlim[2] - xlim[1]), alim[2] - (1.4/100) * 
                  cex * (xlim[2] - xlim[1]), alim[2]), y = c(-1, 
                  -1 + (height/150) * cex * efac[1], -1 - (height/150) * 
                    cex * efac[1], -1), col = col[2], border = col[2], 
                  ...)
            }
        }
        polygon(x = c(b.ci.lb, b, b.ci.ub, b), y = c(-1, -1 + 
            (height/100) * cex * efac[2], -1, -1 - (height/100) * 
            cex * efac[2]), col = col[1], border = border, ...)
        if (missing(mlab)) 
            mlab <- ifelse((x$method == "FE"), "FE Model", "RE Model")
        text(xlim[1], -1, mlab, pos = 4, cex = cex, ...)
    }
    axis(side = 1, at = at, labels = at.lab, cex.axis = cex.axis, 
        ...)
    if (missing(xlab)) 
        xlab <- .setlab(measure, transf.char, atransf.char, gentype = 1)
    mtext(xlab, side = 1, at = min(at) + (max(at) - min(at))/2, 
        line = par("mgp")[1] - 0.5, cex = cex.lab, ...)
    for (i in seq_len(k)) {
        if (is.na(yi[i]) || is.na(vi[i])) 
            next
        if (ci.lb[i] >= alim[2]) {
            polygon(x = c(alim[2], alim[2] - (1.4/100) * cex * 
                (xlim[2] - xlim[1]), alim[2] - (1.4/100) * cex * 
                (xlim[2] - xlim[1]), alim[2]), y = c(rows[i], 
                rows[i] + (height/150) * cex * efac[1], rows[i] - 
                  (height/150) * cex * efac[1], rows[i]), col = "black", 
                ...)
            next
        }
        if (ci.ub[i] <= alim[1]) {
            polygon(x = c(alim[1], alim[1] + (1.4/100) * cex * 
                (xlim[2] - xlim[1]), alim[1] + (1.4/100) * cex * 
                (xlim[2] - xlim[1]), alim[1]), y = c(rows[i], 
                rows[i] + (height/150) * cex * efac[1], rows[i] - 
                  (height/150) * cex * efac[1], rows[i]), col = "black", 
                ...)
            next
        }
        segments(max(ci.lb[i], alim[1]), rows[i], min(ci.ub[i], 
            alim[2]), rows[i], lty = lty[1], ...)
        if (ci.lb[i] >= alim[1]) {
            segments(ci.lb[i], rows[i] - (height/150) * cex * 
                efac[1], ci.lb[i], rows[i] + (height/150) * cex * 
                efac[1], ...)
        }
        else {
            polygon(x = c(alim[1], alim[1] + (1.4/100) * cex * 
                (xlim[2] - xlim[1]), alim[1] + (1.4/100) * cex * 
                (xlim[2] - xlim[1]), alim[1]), y = c(rows[i], 
                rows[i] + (height/150) * cex * efac[1], rows[i] - 
                  (height/150) * cex * efac[1], rows[i]), col = "black", 
                ...)
        }
        if (ci.ub[i] <= alim[2]) {
            segments(ci.ub[i], rows[i] - (height/150) * cex * 
                efac[1], ci.ub[i], rows[i] + (height/150) * cex * 
                efac[1], ...)
        }
        else {
            polygon(x = c(alim[2], alim[2] - (1.4/100) * cex * 
                (xlim[2] - xlim[1]), alim[2] - (1.4/100) * cex * 
                (xlim[2] - xlim[1]), alim[2]), y = c(rows[i], 
                rows[i] + (height/150) * cex * efac[1], rows[i] - 
                  (height/150) * cex * efac[1], rows[i]), col = "black", 
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
                if (addfit && x$int.only) {
                  annotext <- round(cbind(sapply(c(yi, b), atransf), 
                    sapply(c(ci.lb, b.ci.lb), atransf), sapply(c(ci.ub, 
                      b.ci.ub), atransf)), digits[1])
                }
                else {
                  annotext <- round(cbind(sapply(yi, atransf), 
                    sapply(ci.lb, atransf), sapply(ci.ub, atransf)), 
                    digits[1])
                }
            }
            else {
                if (addfit && x$int.only) {
                  annotext <- round(cbind(sapply(c(yi, b), atransf, 
                    targs), sapply(c(ci.lb, b.ci.lb), atransf, 
                    targs), sapply(c(ci.ub, b.ci.ub), atransf, 
                    targs)), digits[1])
                }
                else {
                  annotext <- round(cbind(sapply(yi, atransf, 
                    targs), sapply(ci.lb, atransf, targs), sapply(ci.ub, 
                    atransf, targs)), digits[1])
                }
            }
            rev.order <- ifelse(annotext[, 3] < annotext[, 2], 
                TRUE, FALSE)
            rev.order[is.na(rev.order)] <- FALSE
            annotext[rev.order, 2:3] <- annotext[rev.order, 3:2]
        }
        else {
            if (addfit && x$int.only) {
                annotext <- round(cbind(c(yi, b), c(ci.lb, b.ci.lb), 
                  c(ci.ub, b.ci.ub)), digits[1])
            }
            else {
                annotext <- round(cbind(yi, ci.lb, ci.ub), digits[1])
            }
        }
        if (showweights) {
            if (addfit && x$int.only) {
                annotext <- cbind(round(c(weights, 100), digits[1]), 
                  annotext)
            }
            else {
                annotext <- cbind(round(weights, digits[1]), 
                  annotext)
            }
            annotext <- matrix(apply(annotext, 2, format, nsmall = digits[1]), 
                ncol = 4)
            annotext <- cbind(annotext[, 1], "%    ", annotext[, 
                2], " [ ", annotext[, 3], " , ", annotext[, 4], 
                " ]")
        }
        else {
            annotext <- matrix(apply(annotext, 2, format, nsmall = digits[1]), 
                ncol = 3)
            annotext <- cbind(annotext[, 1], " [ ", annotext[, 
                2], " , ", annotext[, 3], " ]")
        }
        annotext <- apply(annotext, 1, paste, collapse = "")
        if (addfit && x$int.only) {
            text(x = xlim[2], c(rows, -1), labels = annotext, 
                pos = 2, cex = cex, ...)
        }
        else {
            text(x = xlim[2], rows, labels = annotext, pos = 2, 
                cex = cex, ...)
        }
    }
    for (i in seq_len(k)) {
        if (is.na(yi[i])) 
            next
        if (yi[i] >= alim[1] && yi[i] <= alim[2]) 
            points(yi[i], rows[i], pch = pch[i], cex = cex * 
                psize[i], ...)
    }
    if (x$int.only && addfit) 
        abline(h = 0, lty = lty[3], ...)
    res <- list(xlim = par("usr")[1:2], alim = alim, at = at, 
        ylim = ylim, rows = rows, cex = cex, cex.lab = cex.lab, 
        cex.axis = cex.axis)
    invisible(res)
}
