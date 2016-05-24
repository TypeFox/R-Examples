
"plot.predarfima" <- function(x, xlab = NULL, ylab = NULL, main = NULL, ylim = NULL, numback = 5, 
    xlim = NULL, ...) {
    
    op <- par(no.readonly = TRUE)
    z <- x$z
    n <- length(z)
    if (length(xlab) == 0) 
        xlab <- "Time"
    if (length(ylab) == 0) 
        ylab <- "Mode "
    if (length(main) == 0) 
        main <- paste("Time Series and Predictions of ", x$name, sep = "")
    m <- length(x) - 7
    n.ahead <- length(as.vector(x[[1]]$Forecast))
    seed <- x$seed
    limiting <- x$limiting
    bootpred <- x$bootpred
    B <- x$B
    predint <- x$predint
    
    numSD <- qnorm(1 - (1 - predint)/2, sd = 1)
    
    if (numback > 0) {
        ys <- z[(n - numback + 1):n]
        xs <- (n - numback + 1):n
    } else {
        cat("setting numback to 1\n")
        ys <- z[n]
        xs <- n
    }
    islim <- !is.null(x[[1]]$limitSD)
    
    minn <- if (numback > 0) 
        min(ys) else Inf
    maxx <- if (numback > 0) 
        max(ys) else -Inf
    for (i in 1:m) {
        exactPIU <- x[[i]]$Forecast + numSD * x[[i]]$exactSD
        exactPIL <- x[[i]]$Forecast - numSD * x[[i]]$exactSD
        limitPIU <- if (islim) 
            x[[i]]$Forecast + numSD * x[[i]]$limitSD else NULL
        limitPIL <- if (islim) 
            x[[i]]$Forecast - numSD * x[[i]]$limitSD else NULL
        maxx <- max(maxx, exactPIU, limitPIU, x[[i]]$uppernp)
        minn <- min(minn, exactPIL, limitPIL, x[[i]]$lowernp)
        x[[i]]$exactPIU <- exactPIU
        x[[i]]$exactPIL <- exactPIL
        x[[i]]$limitPIU <- limitPIU
        x[[i]]$limitPIL <- limitPIL
        x[[i]]$B <- x$B
        x[[i]]$predint <- predint
    }
    ranger <- maxx - minn
    if (length(ylim) == 0) 
        ylim <- c(minn - ranger/20, maxx + ranger/20)
    
    leg <- c("Exact prediction")
    ltt <- c(1)
    coll <- c("gray")
    
    ltt <- c(ltt, 2)
    coll <- c(coll, "red")
    leg <- c(leg, paste("Exact", predint * 100, "% PI"))
    if (islim) {
        ltt <- c(ltt, 2)
        coll <- c(coll, "orange")
        leg <- c(leg, paste("Limiting", predint * 100, "% PI"))
    }
    
    if (!is.null(x[[2]]$uppernp)) {
        ltt <- c(ltt, 1, 2)
        coll <- c(coll, "blue", "blue")
        leg <- c(leg, paste("Bootstrap prediction"), paste("Bootstrap", predint * 100, "% PI"))
    }
    nf <- layout(matrix(c(1:(m + 1)), m + 1, 1), widths = rep(1, m + 1), heights = c(0.5, 
        rep(1, m)), FALSE)
    par(mar = c(0, 0, 0, 0))
    plot.default(0, 0, axes = FALSE, type = "n")
    legend("bottom", legend = leg, col = coll, lty = ltt, bty = "n", ncol = ceiling(length(leg)/2))
    text(0, 0.5, main, cex = 2)
    for (i in 1:m) {
        par(mar = c(3, 4.2, 0, 0))
        xx <- x[[i]]
        xx$z <- ys
        xx$xs <- c(xs, xs[length(xs)] + 1:(n.ahead + 1))
        ylab1 <- paste(ylab, i, sep = "")
        plotpredARFIMA(xx, xlab = xlab, ylab = ylab1, ylim = ylim)
    }
    par(op)
}




"plotpredARFIMA" <- function(xy, xlab = NULL, ylab = NULL, ylim = NULL, ...) {
    nn <- length(xy$z)
    plot(x = xy$xs[1:nn], y = xy$z, ylab = ylab, xlab = xlab, xlim = c(xy$xs[1], xy$xs[length(xy$xs)]), 
        ylim = ylim, type = "l", cex.lab = 1.4)
    zlast <- xy$z[nn]
    lines(xy$xs[nn:(length(xy$xs) - 1)], c(zlast, xy$Forecast), col = "gray")
    lines(xy$xs[(nn + 1):(length(xy$xs) - 1)], xy$exactPIU, col = "red", lty = 2)
    lines(xy$xs[(nn + 1):(length(xy$xs) - 1)], xy$exactPIL, col = "red", lty = 2)
    predint <- xy$predint
    
    if (!is.null(xy$limitPIU)) {
        lines(xy$xs[(nn + 1):(length(xy$xs) - 1)], xy$limitPIU, col = "orange", lty = 2)
        lines(xy$xs[(nn + 1):(length(xy$xs) - 1)], xy$limitPIL, col = "orange", lty = 2)
    }
    if (!is.null(xy$uppernp)) {
        lines(xy$xs[nn:(length(xy$xs) - 1)], c(zlast, xy$meanvalnp), col = "blue")
        lines(xy$xs[(nn + 1):(length(xy$xs) - 1)], xy$uppernp, col = "blue", lty = 2)
        lines(xy$xs[(nn + 1):(length(xy$xs) - 1)], xy$lowernp, col = "blue", lty = 2)
    }
    
} 
