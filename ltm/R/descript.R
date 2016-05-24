descript <-
function (data, n.print = 10, chi.squared = TRUE, B = 1000) {
    if (!inherits(data, "matrix") && !inherits(data, "data.frame"))
        stop("'data' must be either a data.frame or a matrix")
    nam <- deparse(substitute(data))
    perc <- if (is.matrix(data)) {
        apply(data, 2, function (x) as.vector(table(x)) / sum(!is.na(x)))
    } else {
        out <- lapply(data, function (x) {
            tab <- table(x)
            nams <- names(tab)
            out <-  tab / sum(!is.na(x))
            names(out) <- nams
            out
        })
        nams <- lapply(out, names)
        nams1 <- nams[[1]]
        if (all(sapply(nams, function (x) isTRUE(all.equal(x, nams1)))))
            out <- matrix(unlist(out), length(nams[[1]]), dimnames = list(nams1, names(out)))
        out
    }
    X <- data.matrix(data)
    n <- nrow(X)
    p <- ncol(X)
    missin <- if (any(is.na(X))) {
        apply(X, 2, function (x) {
            nas <- sum(is.na(x))
            c(Freq = nas, "(%)" = 100 * nas / n)
        })
    } else
        NULL
    dat <- na.exclude(data)    
    attr(dat, "na.action") <- NULL    
    ind <- t(combn(p, 2))
    nind <- nrow(ind)
    pvals <- numeric(nind)
    old <- options(warn = (-1))
    on.exit(options(old))
    if (chi.squared) {
        for (i in 1:nind) {
            chsq <- chisq.test(table(X[, ind[i, 1]], X[, ind[i, 2]]))
            if (any(chsq$expected < 5))
                chsq <- chisq.test(table(X[, ind[i, 1]], X[, ind[i, 2]]), simulate.p.value = TRUE, B = B)
            pvals[i] <- chsq$p.value
        }
        pw.ass <- data.frame(ind, format.pval(pvals, digits = 1))
        names(pw.ass) <- c("Item i", "Item j", "p.value")
        pw.ass <- pw.ass[order(pvals, decreasing = TRUE), ]
        row.names(pw.ass) <- 1:nind
    }
    levs. <- apply(dat, 2, function (x) length(unique(x)))
    levs <- pmax(2, levs.)
    levs <- if (all(levs == 2)) 0:sum(levs - 1) else p:sum(levs)
    totSc <- rowSums(X)
    itms <- rbind(Freq = table( factor(totSc, levels = levs) ))
    out <- list(sample = c(p, n), perc = if (is.matrix(perc)) t(perc) else perc, items = itms, 
                pw.ass = if (chi.squared) pw.ass else NULL, n.print = n.print, name = nam, missin = missin, data = data)
    if (all(levs. == 2)) {
        out$bisCorr <- apply(X, 2, biserial.cor, x = totSc, use = "complete.obs", level = 2)
        res <- numeric(p)
        for (i in 1:p) {
            res[i] <- biserial.cor(rowSums(X[, -i]), X[, i], use = "complete.obs", level = 2)
        }
        out$ExBisCorr <- res
        out$perc <- cbind(out$perc, "logit" = qlogis(out$perc[, 2]))
    }
    alpha <- cbind(c(cronbach.alpha(X, na.rm = TRUE)$alpha, 
                    sapply(1:p, function (i) cronbach.alpha(X[, -i], na.rm = TRUE)$alpha)))
    nams <- if (is.null(cnams <- colnames(X))) paste("Excluding Item", 1:p) else paste("Excluding", cnams)
    dimnames(alpha) <- list(c("All Items", nams), "value")
    out$alpha <- alpha
    class(out) <- "descript"
    out
}
