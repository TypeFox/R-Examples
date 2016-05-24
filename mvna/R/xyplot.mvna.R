xyplot.mvna <- function(x, data = NULL, xlab = "Time", ylab = "Cumulative Hazard",
                        tr.choice = "all", conf.int = TRUE, var.type = c("aalen", "greenwood"),
                        ci.fun = c("log", "linear", "arcsin"), level = 0.95, col = c(1, 1, 1),
                        lty = c(1, 3, 3), ci.type = c(1, 2), ...) {

    if (!inherits(x, "mvna"))
        stop("'x' must be of class 'mvna'")
    if (!(length(ci.type)==2 & ci.type[1] %in% c(1,2) & ci.type[2] %in% c(1,2,3)))
        stop("Argument 'ci.type' is not in the good range. See documentation files")
    if (level < 0 | level > 1)
        stop("Argument 'level' must be between 0 and 1")

    ntrans <- nrow(x$trans)
    names.trans <- paste(x$trans[, 1], x$trans[, 2])
    if (tr.choice[1] == "all") {
        tr.choice <- names.trans
    }
    if (sum(tr.choice %in% names.trans) != length(tr.choice)) 
        stop("'tr.choice' and the name of the possible transitions must match")

    z <- mvna::summary.mvna(x, level = level, var.type = var.type,
                            ci.type = ci.type)

    z <- z[tr.choice]
    zz <- do.call(rbind, z)
    zz$cov <- rapply(mapply(rep, tr.choice, sapply(z, nrow), SIMPLIFY = FALSE), c)
    ## sthg to keep the order of tr.choice
    zz$cov <- factor(zz$cov, levels = tr.choice)
    
    ## Let's do the plots
    if (conf.int) {
        dessin <- xyplot(na + lower + upper ~ time | cov, data = zz, type = "s",
                         col = col, lty = lty,
                         xlab = xlab, ylab = ylab, ...)
    } else {
        dessin <- xyplot(na ~ time | cov, data = zz, type = "s",
                         col = col, lty = lty, xlab = xlab, ylab = ylab, ...)
    }
    
    dessin
}
