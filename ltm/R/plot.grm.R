plot.grm <-
function (x, type = c("ICC", "IIC", "OCCu", "OCCl"), items = NULL, category = NULL, 
                      zrange = c(-3.8, 3.8), z = seq(zrange[1], zrange[2], length = 100), annot, 
                      labels = NULL, legend = FALSE, 
                      cx = "top", cy = NULL, ncol = 1, bty = "n", col = palette(), lty = 1, pch, 
                      xlab, ylab, main, sub = NULL, cex = par("cex"), cex.lab = par("cex.lab"), 
                      cex.main = par("cex.main"), cex.sub = par("cex.sub"), cex.axis = par("cex.axis"), 
                      plot = TRUE, ...) {
    if (!inherits(x, "grm"))
        stop("Use only with 'grm' objects.\n")
    type <- match.arg(type)
    betas <- x$coefficients
    nitems <- length(betas)
    ncatg <- sapply(betas, length)
    itms <- if (!is.null(items)) {
        if (!is.numeric(items) || length(items) > nitems)
            stop("'items' must be a numeric vector of length at most ", nitems)
        if (type == "ICC" && any(items < 1 | items > nitems))
            stop("'items' must contain numbers between 1 and ", nitems, " denoting the items.\n")
        if (type == "IIC" && any(items < 0 | items > nitems))
            stop("'items' must contain numbers between 0 and ", nitems)
        items
    } else
        1:nitems
    ctg <- if (!is.null(category)) {
        if (length(category) > 1)
            stop("'category' must be a number indicating the category.\n")
        if (category < 0 || category > max(ncatg))
            stop(paste("'category' must be a number between 1 and ", max(ncatg), ".\n", sep = ""))
        if (any(ind <- category > ncatg)){
            if (sum(ind) > 1)
                warning("Items ", paste(items[ind], collapse = ", "), " are excluded since they have only ", 
                            paste(ncatg[ind], collapse = ", "), " categories, respectively.\n")
            else
                warning("Item ", items[ind], " is excluded since they have only ", ncatg[ind], " categories.\n")
            itms <- itms[!ind]
        }
        category
    } else
        "all"
    cpr <- switch(type, "ICC" = iprobs(betas, z), "IIC" = infoprobs(betas, z), 
                    "OCCl" = cumprobs(betas, z), "OCCu" = cumprobs(betas, z, FALSE))
    plot.items <- type == "ICC" || (type == "IIC" & (is.null(items) || all(items > 0)))
    plot.info <- !plot.items
    if (missing(main)) {
        Main <- if (type == "ICC") {
            "Item Response Category Characteristic Curves"
        } else if (type == "OCCl" || type == "OCCu") {
            "Item Operation Characteristic Curves"
        } else {
            if (plot.items) "Item Information Curves" else "Test Information Function"
        }
        mis.ind <- TRUE
    } else 
        mis.ind <- FALSE
    if (missing(ylab)) {
        ylab <- if (type == "ICC" || (type == "OCCl" | type == "OCCu")) "Probability" else "Information"
    }
    if (missing(xlab)) {
        xlab <- "Ability"
    }
    if (missing(annot)) {
        annot <- !legend
    }
    col. <- col
    lty. <- lty
    if (type == "ICC" || (type == "OCCl" | type == "OCCu")) {
        if (ctg == "all") {
            if (plot) {
                if (length(itms) > 1) {
                    old.par <- par(ask = TRUE)
                    on.exit(par(old.par))
                }
                one.fig <- prod(par("mfcol")) == 1
                for (i in seq(along = itms)) {
                    ii <- itms[i]
                    if (mis.ind) {
                        main. <- if (one.fig) paste("- Item:", names(cpr)[ii]) else paste("\nItem:", names(cpr)[ii])
                        main <- paste(Main, main.)
                    }
                    plot(range(z), c(0, 1), type = "n", xlab = xlab, ylab = ylab, main = main, sub = sub, cex = cex, 
                         cex.lab = cex.lab, cex.main = cex.main, cex.axis = cex.axis, cex.sub = cex.sub, ...)
                    p <- cpr[[ii]]
                    pos <- round(seq(10, 90, length = ncol(p)))
                    col <- rep(col., length.out = ncatg[ii])
                    lty <- rep(lty., length.out = ncatg[ii])
                    if (!missing(pch)) {
                        pch <- rep(pch, length.out = ncatg[ii])
                        pch.ind <- round(seq(15, 85, length = 4))
                    }
                    for (j in 1:ncol(p)) {
                        lines(z, p[, j], lty = lty[j], col = col[j], ...)
                        if (!missing(pch))
                            points(z[pch.ind], p[pch.ind, j], pch = pch[j], col = col[j], cex = cex, ...)
                        if (annot)
                            text(z[pos[j]], p[pos[j], j], adj = c(0, 1.2), if (missing(labels)) j else labels[j], 
                                 col = col[j], cex = cex, ...)
                    }
                    if (legend) {
                        ncol. <- if (is.null(ncol)) {
                            if (is.factor(x$X[[ii]]) && any(nchar(levels(x$X[[ii]])) > 11)) 1 else ncatg[ii] %/% 2
                        } else
                            ncol
                        lab <- if (missing(labels)) {
                            switch(type,
                                "ICC" = if (is.factor(x$X[[ii]])) levels(x$X[[ii]]) else 1:ncatg[ii],
                                "OCCu" = paste(if (is.factor(x$X[[ii]])) levels(x$X[[ii]])[-1] else 2:ncatg[ii], "or higher"),
                                "OCCl" = paste(if (is.factor(x$X[[ii]])) levels(x$X[[ii]])[-1] else 2:ncatg[ii], "or lower"))
                        } else
                            labels
                        legend(cx, cy, legend = lab, lty = lty, pch = pch, col = col, bty = bty, ncol = ncol., cex = cex, ...)
                    }
                }
                return.value <- list(z = z, pr = cpr[itms])
            } else {
                return(list(z = z, pr = cpr[itms]))
            }
        } else {
            p <- sapply(cpr[itms], function (x, category) x[, category], category = ctg)
            if (plot) {
                one.fig <- prod(par("mfcol")) == 1
                if (mis.ind) {
                    main. <- if (one.fig) paste("- Category:", ctg) else paste("\nCategory:", ctg)
                    main <- paste(Main, main.)
                }
                plot(range(z), c(0, 1), type = "n", xlab = xlab, ylab = ylab, main = main, sub = sub, cex = cex, 
                     cex.lab = cex.lab, cex.main = cex.main, cex.axis = cex.axis, cex.sub = cex.sub, ...)
                pos <- round(seq(10, 90, length = ncol(p)))
                col <- rep(col., length.out = length(itms))
                lty <- rep(lty., length.out = length(itms))
                if (!missing(pch)) {
                    pch <- rep(pch, length.out = length(itms))
                    pch.ind <- round(seq(15, 85, length = 4))
                }                   
                for (j in 1:ncol(p)) {
                    lines(z, p[, j], lty = lty[j], col = col[j], ...)
                    if (!missing(pch))
                        lines(z, p[, j], pch = pch[j], col = col[j], cex = cex, ...)
                    if (annot)
                        text(z[pos[j]], p[pos[j], j], adj = c(0, 1.2), if (missing(labels)) colnames(p)[j] else labels[j], 
                             col = col[j], cex = cex, ...)
                }
                if (legend) {
                    ncol. <- if (is.null(ncol)) {
                        if (any(nchar(colnames(p)) > 11)) 1 else length(items) %/% 2
                    } else 
                        ncol
                    lab <- if (missing(labels)) colnames(p) else labels
                    legend(cx, cy, legend = lab, lty = lty, pch = pch, col = col, bty = bty, ncol = ncol., cex = cex, ...)
                }
                return.value <- cbind(z = z, pr = p)
            } else {
                return(cbind(z = z, pr = p))
            }
        }
    } else {
        p <- cpr[, itms, drop = FALSE]
        if (plot) {
            r <- if (plot.items) range(p) else range(rowSums(cpr))
            if (mis.ind) {
                main <- Main
            }
            plot(range(z), r, type = "n", xlab = xlab, ylab = ylab, main = main, sub = sub, cex = cex, 
                 cex.lab = cex.lab, cex.main = cex.main, cex.axis = cex.axis, cex.sub = cex.sub, ...)
            if (plot.items) {
                col <- rep(col., length.out = length(itms))
                lty <- rep(lty., length.out = length(itms))
                if (!missing(pch)) {
                    pch <- rep(pch, length.out = length(itms))
                    pch.ind <- round(seq(15, 85, length = 4))
                }
                pos <- round(seq(10, 90, length = ncol(p)))
                for (i in seq(along = itms)) {
                    lines(z, p[, i], lty = lty[i], col = col[i], ...)
                    if (!missing(pch))
                        points(z[pch.ind], p[pch.ind, i], pch = pch[i], col = col[i], cex = cex, ...)
                    if (annot)
                        text(z[pos[i]], p[pos[i], i], adj = c(0, 1.2), if (missing(labels)) i else labels[i], 
                             col = col[i], cex = cex, ...)
                }
                if (legend) {
                    ncol. <- if (is.null(ncol)) ncol. <- if (nitems > 8) 2 else 1 else ncol
                    legend(cx, cy, legend = if (missing(labels)) colnames(cpr)[itms] else labels, lty = lty, pch = pch, 
                           col = col, bty = bty, ncol = ncol., cex = cex, ...)
                }
            } else {
                col <- col.[1]
                lty <- lty.[1]
                p <- rowSums(cpr)
                lines(z, p, lty = lty, col = col, ...)
                if (!missing(pch))
                    points(z, p, pch = pch, col = col, cex = cex, ...)
                if (legend)
                    legend(cx, cy, legend = "Information", lty = lty, col = col, pch = pch, bty = bty, ncol = 1, ...)
            }
            return.value <- if (plot.items) cbind(z = z, item.info = p) else cbind(z = z, test.info = p)
        } else {
            return(if (plot.items) cbind(z = z, item.info = p) else cbind(z = z, test.info = rowSums(cpr)))
        }
    }
    invisible(return.value)
}
