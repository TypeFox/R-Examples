plot.tpm <-
function (x, type = c("ICC", "IIC"), items = NULL, zrange = c(-3.8, 3.8),
                      z = seq(zrange[1], zrange[2], length = 100), annot, labels = NULL, 
                      legend = FALSE, cx = "topleft", cy = NULL, ncol = 1, bty = "n", col = palette(), lty = 1, 
                      pch, xlab, ylab, main, sub = NULL, cex = par("cex"), cex.lab = par("cex.lab"), 
                      cex.main = par("cex.main"), cex.sub = par("cex.sub"), cex.axis = par("cex.axis"), 
                      plot = TRUE, ...) {
    if (!inherits(x, "tpm"))
        stop("Use only with 'tpm' objects.\n")
    type <- match.arg(type)
    thetas <- x$coefficients
    cs <- plogis(thetas[, 1]) * x$max.guessing
    betas <- thetas[, 2:3]
    p <- nrow(betas)
    itms <- if (!is.null(items)) {
        if (!is.numeric(items) || length(items) > p)
            stop("'items' must be a numeric vector of length at most ", p)
        if (type == "ICC" && any(items < 1 | items > p))
            stop("'items' must contain numbers between 1 and ", p, " denoting the items.\n")
        if (type == "IIC" && any(items < 0 | items > p))
            stop("'items' must contain numbers between 0 and ", p)
        items
    } else
        1:p
    Z <- cbind(1, z)
    pr <- if (type == "ICC") {
        cs <- matrix(cs, length(z), p, TRUE)
        cs + (1 - cs) * probs(Z %*% t(betas))
    } else {
        pi. <- plogis(Z %*% t(betas))
        cs <- matrix(cs, length(z), p, TRUE)
        pi <- cs + (1 - cs) * pi.
        pqr <- pi * (1 - pi) * (pi. / pi)^2
        t(t(pqr) * betas[, 2]^2)
    }
    plot.items <- type == "ICC" || (type == "IIC" & (is.null(items) || all(items > 0)))
    plot.info <- !plot.items
    if (plot) {
        col <- if (plot.items) rep(col, length.out = length(itms)) else col[1]
        lty <- if (plot.items) rep(lty, length.out = length(itms)) else lty[1]
        if(!missing(pch)) {
            pch <- if (plot.items) rep(pch, length.out = length(itms)) else pch[1]
            pch.ind <- round(seq(15, 85, length = 4))
        }
        if (missing(main)) {
            main <- if (type == "ICC") "Item Characteristic Curves" else { 
                if (plot.items) "Item Information Curves" else "Test Information Function"
            }
        }
        if (missing(xlab)) {
            xlab <- "Ability"
        }
        if (missing(ylab)) {
            ylab <- if (type == "ICC") "Probability" else "Information"
        }
        r <- if (type == "ICC") c(0, 1) else { if (plot.info) range(rowSums(pr)) else range(pr[, itms]) }
        plot(range(z), r, type = "n", xlab = xlab, ylab = ylab, main = main, sub = sub, cex = cex, cex.lab = cex.lab, 
             cex.main = cex.main, cex.axis = cex.axis, cex.sub = cex.sub, ...)
        if (missing(annot)) {
            annot <- !legend
        }
        if (legend) {
            legnd <- if (is.null(labels)) {
                if (plot.info) "Information" else rownames(betas)[itms]
            } else {
                if (length(labels) < length(itms))
                    warning("the length of 'labels' is smaller than the length of 'items'.\n")
                labels
            }
            legend(cx, cy, legend = legnd, col = col, lty = lty, bty = bty, ncol = ncol, cex = cex, pch = pch, ...) 
        }
        if (annot) {
            pos <- round(seq(10, 90, length = length(itms)))
            nams <- if (is.null(labels)) { 
                nms <- if (rownames(betas)[1] == "Item 1") 1:p else rownames(betas)
                nms[itms]
            } else {
                if (length(labels) < length(itms))
                    warning("the length of 'labels' is smaller than the length of 'items'.\n")
                labels
            }
        }
        if (plot.items) {
            for (it in seq(along = itms)) {
                lines(z, pr[, itms[it]], lty = lty[it], col = col[it], ...)
                if (!missing(pch))
                    points(z[pch.ind], pr[pch.ind, itms[it]], pch = pch[it], col = col[it], cex = cex, ...)
                if (annot)
                    text(z[pos[it]], pr[pos[it], itms[it]], labels = nams[it], adj = c(0, 2), col = col[it], cex = cex, ...)
            }
        }
        if (plot.info)
            lines(z, rowSums(pr), lty = lty, col = col, ...)
        invisible(if (plot.items) cbind(z = z, pr[, itms]) else cbind(z = z, info = rowSums(pr)))
    } else {
        if (plot.items) cbind(z = z, pr[, itms]) else cbind(z = z, info = rowSums(pr))
    }
}
