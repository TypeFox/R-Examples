lines.mvna <- function(x, tr.choice, col = 1, lty,
                       conf.int = FALSE, level = 0.95,
                       var.type = c("aalen", "greenwood"),
                       ci.fun = c("log", "linear", "arcsin"),
                       ci.col = col, ci.lty = 3, ...) {

    if (!inherits(x, "mvna")) {
        stop("'x' must be of class 'mvna'")
    }
        ref <- paste(x$trans[, 1], x$trans[, 2])
    if (missing(tr.choice)) {
        tr.choice <- ref
    }
    if (sum(tr.choice %in% ref) != length(tr.choice)) {
        stop("Names of the possible transitions and 'tr.choice' must match")
    }

    object <- mvna::summary.mvna(x, level = level, var.type = var.type,
                                 ci.fun = ci.fun)[tr.choice]
    lt <- length(object)

    if (missing(lty)) {
        lty <- seq_len(lt)
    }
    else if (length(lty) < lt) {
        lty <- lty * rep(1, lt)
    }
    if (length(col) < lt)
        col <- col * rep(1, lt)

        for (i in seq_len(lt)) {
        lines(object[[i]]$time, object[[i]]$na, type = "s",
              col = col[i], lty = lty[i], ...)
    }

    if (conf.int) {
        if (length(ci.col) < lt)
            ci.col <- ci.col * rep(1, lt)
        if (length(ci.lty) < lt)
            ci.lty <- ci.lty * rep(1, lt)
        for (i in seq_len(lt)) {
            lines(object[[i]]$time, object[[i]]$lower, type = "s",
                  col = ci.col[i], lty = ci.lty[i], ...)
            lines(object[[i]]$time, object[[i]]$upper, type = "s",
                  col = ci.col[i], lty = ci.lty[i], ...)
        }
    }
    invisible()
}
