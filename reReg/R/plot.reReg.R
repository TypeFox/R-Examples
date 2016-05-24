plot.reReg <- function(x, se = FALSE, B = 200, breaks = 1000, ...) {
    if (!is.reReg(x))
        stop("Response must be a reReg xect")
    ly <- hy <- lyU <- lyL <- hyU <- hyL <- NULL
    alpha <- x$alpha
    beta <- x$beta
    id <- x$df$id
    cluster <- unlist(lapply(split(id, id), function(z) 1:length(z)))
    clsz <- unlist(lapply(split(id, id), length))
    Y <- rep(x$df$Time[cumsum(clsz)], clsz)
    T <- x$df$Time
    X <- as.matrix(x$df[,-c(1:3)])
    delta <- x$df$event
    t <- seq(0, max(Y), length.out = breaks)
    Ya <- log(Y) + X %*% alpha
    Ta <- log(T) + X %*% alpha
    mt <- unlist(lapply(split(cluster, id), length)) - 1
    lambda <- npMLE(Ya[which(cluster == 1)], Ta, Ya)
    ly <- npMLE(t, exp(Ta), exp(Ya))
    zHat <- as.numeric(mt * max(ly) / lambda)
    ly <- ly / max(ly)
    win.ly <- max(ly)
    if (se) {
        E <- matrix(rexp(length(t) * B), nrow = length(t))
        lytmp <- apply(E, 2, function(x) npMLE(t, exp(Ta), exp(Ya), x))
        lytmp <- apply(lytmp, 2, function(z) z / max(z))
        lyU <- apply(lytmp, 1, function(z) quantile(z, 0.975))
        lyL <- apply(lytmp, 1, function(z) quantile(z, 0.025))
        win.ly <- max(lyU)
    }
    zHat <- ifelse(zHat %in% c("Inf", "NA", "NaN"), 0, zHat)
    Yb <- log(Y) + X %*% beta
    Yb <- Yb[which(cluster == 1)]
    hy <- sapply(t, function(z) baseHaz(z, exp(Yb), zHat, delta[cluster == 1]))
    win.hy <- max(hy)
    if (se) {
        E <- matrix(rexp(sum(cluster == 1) * B), nrow = sum(cluster == 1))
        hytmp <- apply(E, 2, function(z)
            sapply(t, function(y)
                baseHaz(y, exp(Yb), zHat, delta[which(cluster == 1)], z)))
        hyU <- apply(hytmp, 1, function(z) quantile(z, 0.975))
        hyL <- apply(hytmp, 1, function(z) quantile(z, 0.025))
        win.hy <- max(hyU)
    }
    ## plot(t, ly, type = "l", xlab = "", ylab = "", xaxt = "n", yaxt = "n", ylim = c(0, win.ly),
    ##      main = "Baseline Cumulative Rate Function")
    plot(t, ly, type = "l",  xlab = "", ylab = "", ylim = c(0, win.ly), main = "Baseline Cumulative Rate Function")
    if (se) {
        lines(t, lyU, col = 2)
        lines(t, lyL, col = 2)
    }
    title(ylab = expression(hat(Lambda)[0](t)), xlab = "time", line = 2.2)
    ##     title(xlab = "time", line = 1.6, cex.lab = 0.8)
    ## axis(1, at = seq(0, round(max(t), 1), length.out = 11), 
    ##      cex.axis = 0.8, tck = -0.015, mgp = c(3, .3, 0))
    ## axis(2, at = seq(0, round(win.ly, 2), length.out = 11),
    ##      las = 2, cex.axis = 0.8, tck = -0.015, mgp = c(3, .4, 0))
    op <- par(ask=TRUE)
    ## plot(t, hy, type = "l", xlab = "", ylab = "", xaxt = "n", yaxt = "n", ylim = c(0, win.hy),
    ##      main = "Baseline Cumulative Hazard Function")
    plot(t, hy, type = "l",  xlab = "", ylab = "", ylim = c(0, win.hy), main = "Baseline Cumulative Hazard Function")
    if (se) {
        lines(t, hyU, col = 2, lty = 2)
        lines(t, hyL, col = 2, lty = 2)
    }
    title(ylab = expression(hat(H)[0](t)), xlab = "time", line = 2.2)
    ## title(xlab = "time", line = 1.6, cex.lab = 0.8)
    ## axis(1, at = seq(0, round(max(t), 1), length.out = 11), 
    ##      cex.axis = 0.8, tck = -0.015, mgp = c(3, .3, 0))
    ## axis(2, at = seq(0, round(win.hy, 2), length.out = 11),
    ##      las = 2, cex.axis = 0.8, tck = -0.015, mgp = c(3, .4, 0))
    par(op)
    out <- list(x = t, ly = ly, lyU = lyU, lyL = lyL, hy = hy, hyU = hyU, hyL = hyL)
    invisible(out)
}
