lines.etm <- function(x, tr.choice, col = 1, lty,
                      conf.int = FALSE, level = 0.95, ci.fun = "linear",
                      ci.col = col, ci.lty = 3, ...) {

    if (!inherits(x, "etm")) {
        stop("'x' must be of class 'etm'")
    }
    
    ufrom <- unique(x$trans$from)
    uto <- unique(x$trans$to)
    absorb <- setdiff(uto, ufrom)
    nam1 <- dimnames(x$est)[[1]]
    nam2 <- dimnames(x$est)[[2]]
    pos <- c(paste(nam1[!(nam1 %in% as.character(absorb))],
                   nam2[!(nam2 %in% as.character(absorb))]),
             paste(x$trans$from, x$trans$to))
    if (missing(tr.choice)) tr.choice <- pos
    
    ref <- sapply(1:length(x$state.names), function(i) {
        paste(x$state.names, x$state.names[i])
    })
    ref <- matrix(ref)
    if (sum(tr.choice %in% ref == FALSE) > 0)
        stop("Argument 'tr.choice' and possible transitions must match")

    temp <- ci.transfo(x, tr.choice, level, ci.fun)
    
    lt <- length(temp)

    if (missing(lty)) {
        lty <- seq_len(lt)
    }
    else if (length(lty) < lt) {
        lty <- lty * rep(1, lt)
    }
    if (length(col) < lt)
        col <- col * rep(1, lt)

        for (i in seq_len(lt)) {
        lines(temp[[i]]$time, temp[[i]]$P, type = "s",
              col = col[i], lty = lty[i], ...)
    }

    if (conf.int && !is.null(x$cov)) {
        if (length(ci.col) < lt)
            ci.col <- ci.col * rep(1, lt)
        if (length(ci.lty) < lt)
            ci.lty <- ci.lty * rep(1, lt)
        for (i in seq_len(lt)) {
            lines(temp[[i]]$time, temp[[i]]$lower, type = "s",
                  col = ci.col[i], lty = ci.lty[i], ...)
            lines(temp[[i]]$time, temp[[i]]$upper, type = "s",
                  col = ci.col[i], lty = ci.lty[i], ...)
        }
    }
    
    invisible()
}
