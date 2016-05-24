xyplot.etm <- function(x, data = NULL, tr.choice, col = c(1, 1, 1), lty = c(1, 3, 3),
                       xlab="Time", ylab = "Transition probability",
                       conf.int = TRUE, ci.fun = "linear", level = 0.95, ...) {
    
    if (!inherits(x, "etm"))
        stop("Argument 'x' must be of class 'etm'")
    
    ref <- sapply(1:length(x$state.names), function(i) {
        paste(x$state.names, x$state.names[i])
    })
    ref <- matrix(ref)
    
    if (missing(tr.choice)) {
        ufrom <- unique(x$trans$from)
        uto <- unique(x$trans$to)
        absorb <- setdiff(uto, ufrom)
        nam1 <- dimnames(x$est)[[1]]
        nam2 <- dimnames(x$est)[[2]]
        pos <- c(paste(nam1[!(nam1 %in% as.character(absorb))],
                       nam2[!(nam2 %in% as.character(absorb))]),
                 paste(x$trans$from, x$trans$to))
        tr.choice <- pos
    }
    
    if (sum(tr.choice %in% ref == FALSE) > 0)
        stop("Argument 'tr.choice' and possible transitions must match")
    
    temp <- ci.transfo(x, tr.choice, level, ci.fun)
    
    for (i in seq_along(temp)) {
        temp[[i]]$cov <- names(temp)[i]
    }
    temp <- do.call(rbind, temp)
    temp$cov <- factor(temp$cov, levels = tr.choice)
    
    if (conf.int) {
        aa <- xyplot(temp$P + temp$lower + temp$upper ~ temp$time | temp$cov,
                     type = "s", col = col, lty = lty, xlab = xlab, ylab = ylab, ...)
    }
    else {
        aa <- xyplot(temp$P ~ temp$time | temp$cov, type = "s",
                     col = col, lty = lty, xlab = xlab, ylab = ylab, ...)
    }
    
    aa
}
