plot.ltm <-
function (x, type = c("ICC", "IIC", "loadings"), items = NULL, zrange = c(-3.8, 3.8),
                      z = seq(zrange[1], zrange[2], length = 100), annot, 
                      labels = NULL, legend = FALSE, cx = "topleft", cy = NULL, ncol = 1, bty = "n", 
                      col = palette(), lty = 1, pch, xlab, ylab, zlab, main, sub = NULL, cex = par("cex"), 
                      cex.lab = par("cex.lab"), cex.main = par("cex.main"), cex.sub = par("cex.sub"), 
                      cex.axis = par("cex.axis"), plot = TRUE, ...) {
    if (!inherits(x, "ltm"))
        stop("Use only with 'ltm' objects.\n")
    type <- match.arg(type)
    if (type == "IIC" && any(x$ltst$factors == 2, x$ltst$inter, x$ltst$quad.z1, x$ltst$quad.z2))
        stop("Item Information Curves are currently plotted only for the one-factor model.\n")
    if (type == "loadings" && x$ltst$factors == 1)
        stop("type = 'loadings' is only valid for the linear two-factor model.\n")
    betas <- x$coefficients
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
    if (x$ltst$factors == 1){
        Z <- if (x$ltst$quad.z1) cbind(1, z, z * z) else cbind(1, z)
        pr <- if (type == "ICC") {
            plogis(Z %*% t(betas))
        } else {
            pr <- plogis(Z %*% t(betas))
            pqr <- pr * (1 - pr)
            t(t(pqr) * betas[, 2]^2)
        }
        plot.items <- type == "ICC" || (type == "IIC" & (is.null(items) || all(items > 0)))
        plot.info <- !plot.items
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
        if (plot) {
            plot(range(z), r, type = "n", xlab = xlab, ylab = ylab, main = main, sub = sub, cex = cex, 
                 cex.lab = cex.lab, cex.main = cex.main, cex.axis = cex.axis, cex.sub = cex.sub, ...)
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
                        text(z[pos[it]], pr[pos[it], itms[it]], labels = nams[it], adj = c(0, 2), col = col[it], 
                             cex = cex, ...)
                }
            }
            if (plot.info)
                lines(z, rowSums(pr), lty = lty, col = col, ...)
            return.value <- if (plot.items) cbind(z = z, pr[, itms]) else cbind(z = z, info = rowSums(pr))
        } else {
            return(if (plot.items) cbind(z = z, pr[, itms]) else cbind(z = z, info = rowSums(pr)))
        }
    }
    if (x$ltst$factors == 2) {
        nams <- rownames(x$coef)
        if (type == "loadings") {
            if (any(unlist(x$ltst[2:4])))
                stop("the plot of standardized loadings is produced only for the linear two-factor model.\n")
            cof <- coef(x, TRUE)
            z1 <- cof[itms, 3]
            z2 <- cof[itms, 5]
            if (plot) {
                if (missing(xlab))
                    xlab <- "Factor 1"
                if (missing(ylab))
                    ylab <- "Factor 2"
                if (missing(main))
                    main <- "Standardized Loadings"
                plot(z1, z2, type = "n", xlab = xlab, ylab = ylab, main = main, sub = sub,
                     xlim = c(min(z1, -0.1), max(z1, 0.1)), ylim = c(min(z2, -0.1), max(z2, 0.1)),
                     cex = cex, cex.lab = cex.lab, cex.main = cex.main, cex.axis = cex.axis, cex.sub = cex.sub, ...)
                abline(h = 0, v = 0, lty = 2)
                text(z1, z2, labels = if (is.null(labels)) nams[itms] else labels, cex = cex, ...)
                return.value <- cbind("Factor 1" = z1, "Factor 2" = z2)
            } else {
                return(cbind("Factor 1" = z1, "Factor 2" = z2))
            }
        } else {
            z1 <- z2 <- z
            f <- function (z, betas, strct) {
                Z <- cbind(1, z[1], z[2])
                colnames(Z) <- c("(Intercept)", "z1", "z2")
                if (strct$inter)
                    Z <- cbind(Z, "z1:z2" = z[1] * z[2])
                if (strct$quad.z1)
                    Z <- cbind(Z, "I(z1^2)" = z[1] * z[1])
                if (strct$quad.z2)
                    Z <- cbind(Z, "I(z2^2)" = z[2] * z[2])
                Z <- Z[, match(names(betas), colnames(Z)), drop = FALSE]
                pr <- plogis(Z %*% betas)
            }
            if (plot) {
                if (missing(xlab))
                    xlab <- "Factor 1"
                if (missing(ylab))
                    ylab <- "Factor 2"
                if (missing(zlab))
                    zlab <- "Probability"
                if (missing(main))
                    main <- "Item Characteristic Surfaces"
                col <- if (missing(col)) rep("white", length.out = length(itms)) else rep(col, length.out = length(itms))
                old.par <- par(ask = TRUE)
                on.exit(par(old.par))
            }
            grid. <- as.matrix(expand.grid(z1, z2))
            dimnames(grid.) <- NULL
            z <- vector("list", length(itms))
            for (it in seq(along = itms)) {
                item <- itms[it]
                z[[it]] <- apply(grid., 1, f, betas = betas[item, ], strct = x$ltst)
                dim(z[[it]]) <- c(length(z1), length(z1))
                if (plot) {
                    persp(z1, z2, z[[it]], cex = cex, xlab = list(xlab, cex = cex.lab), 
                          ylab = list(ylab, cex = cex.lab), zlab = list(zlab, cex = cex.lab), 
                          main = list(main, cex = cex.main), sub = list(if (is.null(labels)) nams[item] else labels[it], 
                          cex = cex.sub), col = col[it], ...)
                }
            }
            return.value <- list(Factor1 = z1, Factor2 = z2, Prob = z)
            if (!plot)
                return(return.value)
        }
    }
    invisible(return.value)
}
