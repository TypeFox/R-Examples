kma.show.results <-
function (Result, lwd.functions = 1, lwd.centers = 3) 
{
    if (length(Result) == 0) {
        stop("First parameter (Result) must not be NULL")
    }
    lwd.centers <- lwd.centers
    lwd.functions <- lwd.functions
    iterations <- Result$iterations
    x <- as.matrix(Result$x)
    y0 <- Result$y0
    y1 <- Result$y1
    n.clust <- Result$n.clust
    warping.method <- Result$warping.method
    similarity.method <- Result$similarity.method
    center.method <- Result$center.method
    x.center.orig <- Result$x.center.orig
    y0.center.orig <- Result$y0.center.orig
    y1.center.orig <- Result$y1.center.orig
    similarity.orig <- Result$similarity.orig
    x.final <- Result$x.final
    n.clust.final <- Result$n.clust.final
    x.centers.final <- Result$x.centers.final
    y1.centers.final <- Result$y1.centers.final
    y0.centers.final <- Result$y0.centers.final
    labels <- Result$labels
    similarity.final <- Result$similarity.final
    dilation.list <- Result$dilation.list
    shift.list <- Result$shift.list
    dilation <- Result$dilation
    shift <- Result$shift
    if (length(dim(y1)) != 0) {
        n.camp <- dim(y1)[2]
        n.obs <- dim(y1)[1]
        if (length(dim(y1)) == 3) {
            r <- dim(y1)[3]
        }
        else {
            r <- 1
        }
    }
    if (length(dim(y0)) != 0) {
        n.camp <- dim(y0)[2]
        n.obs <- dim(y0)[1]
        if (length(dim(y0)) == 3) {
            r <- dim(y0)[3]
        }
        else {
            r <- 1
        }
    }
    sim.final <- mean(similarity.final)
    labels.unique <- sort(unique(labels))
    myrainbow <- c("red", "blue", "green3", "orange", "grey", 
        "yellow")
    myrainbow.dark <- c("darkred", "darkblue", "darkgreen", "darkorange", 
        "black", "brown")
    colori.dopo <- rainbow(length(labels.unique))
    myrainbow <- c(myrainbow, colori.dopo)
    myrainbow <- myrainbow[1:n.clust.final]
    myrainbow.dark <- c(myrainbow.dark, colori.dopo)
    myrainbow.dark <- myrainbow.dark[1:n.clust.final]
    colours.random <- rainbow(n.obs)
    colours.bygroup <- rep(0, n.obs)
    for (k in labels.unique) {
        colours.bygroup[which(labels == k)] <- myrainbow[k]
    }
    colours.bygroup.dark <- rep(0, n.obs)
    for (k in labels.unique) {
        colours.bygroup.dark[which(labels == k)] <- myrainbow.dark[k]
    }
    colours.templates.iter1 <- myrainbow
    colours.templates.last <- myrainbow
    colours.warping <- rep(0, n.obs)
    for (k in labels.unique) {
        colours.warping[which(labels == k)] <- myrainbow[k]
        colours.templates.iter1[k] <- myrainbow[k]
        colours.templates.last[k] <- colours.bygroup[which(labels == 
            k)[1]]
    }
    if (length(y0) != 0) {
        if (r == 1) {
            dev.new()
            matplot(t(x), t(y0), type = "l", lwd = lwd.functions, 
                col = colours.random, xlab = "x", ylab = "y")
            title(main = "Original Data")
            for (k in 1:dim(y0.center.orig)[1]) {
                lines(x.center.orig, y0.center.orig[k, ], lwd = lwd.centers, 
                  col = colours.templates.iter1[k])
            }
        }
        else {
            for (l in 1:r) {
                dev.new()
                matplot(t(x), t(y0[, , l]), type = "l", lwd = lwd.functions, 
                  col = colours.random, xlab = "x", ylab = "y")
                title(main = paste(c("Original Data, dimension ", 
                  l), collapse = ""))
                lines(x.center.orig, y0.center.orig[[1]][l, ], 
                  lwd = lwd.centers, col = colours.templates.iter1[l])
            }
        }
    }
    if (length(y0) != 0) {
        if (r == 1) {
            dev.new()
            tex <- paste(c("k = ", n.clust), collapse = "")
            matplot(t(x.final), t(y0), type = "l", lwd = lwd.functions, 
                col = colours.bygroup.dark, xlab = "x", ylab = tex)
            title2 <- c("Registration: ", warping.method)
            title2 <- paste(title2, collapse = "")
            title22 <- c("Aligned Data")
            title2def <- c(title2, title22)
            title(main = title2def)
            for (k in labels.unique) {
                lines(x.centers.final, y0.centers.final[k, ], 
                  lwd = lwd.centers, col = colours.templates.last[k])
            }
            text <- rep(0, length(labels.unique))
            for (i in 1:length(labels.unique)) {
                tt <- c("Cluster ", i)
                tt <- paste(tt, collapse = "")
                text[i] <- tt
            }
            lty <- rep(1, length(text))
            legend("topleft", legend = text, col = colours.templates.last, 
                lty = lty, cex = 0.6)
        }
        else {
            for (l in 1:r) {
                dev.new()
                tex <- paste(c("k = ", n.clust), collapse = "")
                matplot(t(x.final), t(y0[, , l]), type = "l", 
                  lwd = lwd.functions, col = colours.bygroup.dark, 
                  xlab = "x", ylab = tex)
                title2 <- c("Registration: ", warping.method)
                title2 <- paste(title2, collapse = "")
                title22 <- paste(c("Aligned Data; dimension ", 
                  l), collapse = "")
                title2def <- c(title2, title22)
                title(main = title2def)
                for (k in labels.unique) {
                  lines(x.centers.final, y0.centers.final[[k]][l, 
                    ], lwd = lwd.centers, col = colours.templates.last[k])
                }
                text <- rep(0, length(labels.unique))
                for (i in 1:length(labels.unique)) {
                  tt <- c("Cluster ", i)
                  tt <- paste(tt, collapse = "")
                  text[i] <- tt
                }
                lty <- rep(1, length(text))
                legend("topleft", legend = text, col = colours.templates.last, 
                  lty = lty, cex = 0.6)
            }
        }
    }
    if (length(y1) != 0 && (similarity.method == "d1.pearson" || 
        similarity.method == "d1.L2" || similarity.method == 
        "d1.L2.centered")) {
        if (r == 1) {
            dev.new()
            matplot(t(x), t(y1), type = "l", lwd = lwd.functions, 
                col = colours.random, xlab = "x", ylab = "y")
            title(main = "Original Derivative")
            for (k in 1:dim(y1.center.orig)[1]) {
                lines(x.center.orig, y1.center.orig[k, ], lwd = lwd.centers, 
                  col = colours.templates.iter1[k])
            }
        }
        if (r != 1) {
            for (l in 1:r) {
                dev.new()
                matplot(t(x), t(y1[, , l]), type = "l", lwd = lwd.functions, 
                  col = colours.random, xlab = "x", ylab = "y")
                title(main = paste(c("Original Derivative; dimension ", 
                  l), collapse = ""))
                for (k in 1:dim(y1.center.orig)[1]) {
                  lines(x.center.orig, y1.center.orig[[k]][l, 
                    ], lwd = lwd.centers, col = colours.templates.iter1[k])
                }
            }
        }
    }
    if (length(y1) != 0 && (similarity.method == "d1.pearson" || 
        similarity.method == "d1.L2" || similarity.method == 
        "d1.L2.centered")) {
        if (r == 1) {
            dev.new()
            matplot(t(x.final), t(y1), type = "l", lwd = lwd.functions, 
                col = colours.bygroup.dark, xlab = "x", ylab = "y")
            title2 <- c("Registration: ", Result$warping.method)
            title2 <- paste(title2, collapse = "")
            title22 <- c("Aligned Derivative")
            title2def <- c(title2, title22)
            title(main = title2def)
            for (k in labels.unique) {
                lines(x.centers.final, y1.centers.final[k, ], 
                  lwd = lwd.centers, col = colours.templates.last[k])
            }
            text <- rep(0, length(labels.unique))
            for (i in 1:length(labels.unique)) {
                tt <- c("Cluster ", i)
                tt <- paste(tt, collapse = "")
                text[i] <- tt
            }
            lty <- rep(1, length(text))
            legend("topleft", legend = text, col = colours.templates.last, 
                lty = lty, cex = 0.6)
        }
        else {
            for (l in 1:r) {
                dev.new()
                matplot(t(x.final), t(y1[, , l]), type = "l", 
                  lwd = lwd.functions, col = colours.bygroup.dark, 
                  xlab = "x", ylab = "y")
                title2 <- c("Registration: ", Result$warping.method)
                title2 <- paste(title2, collapse = "")
                title22 <- paste(c("Aligned Derivative", l), 
                  collapse = "")
                title2def <- c(title2, title22)
                title(main = title2def)
                for (k in labels.unique) {
                  lines(x.centers.final, y1.centers.final[[k]][l, 
                    ], lwd = lwd.centers, col = colours.templates.last[k])
                }
                text <- rep(0, length(labels.unique))
                for (i in 1:length(labels.unique)) {
                  tt <- c("Cluster ", i)
                  tt <- paste(tt, collapse = "")
                  text[i] <- tt
                }
                lty <- rep(1, length(text))
                legend("topleft", legend = text, col = colours.templates.last, 
                  lty = lty, cex = 0.6)
            }
        }
    }
    dev.new()
    plot(t(x[1, ]), t(x.final[1, ]), xlim = c(min(x, na.rm = TRUE), 
        max(x, na.rm = TRUE)), ylim = c(min(x, na.rm = TRUE), 
        max(x, na.rm = TRUE)), type = "l", lwd = 1, col = colours.warping[1], 
        xlab = "x", ylab = "y", asp = 1)
    for (i in 2:n.obs) {
        lines(t(x[i, ]), t(x.final[i, ]), type = "l", lwd = 1, 
            col = colours.warping[i], xlab = "x", ylab = "y", 
            asp = 1)
    }
    title3 <- c("Registration: ", Result$warping.method)
    title3 <- paste(title3, collapse = "")
    title33 <- c("Warping Functions")
    title3def <- c(title3, title33)
    title(main = title3def)
    abline(v = min(x))
    abline(v = max(x))
    text <- rep(0, length(labels.unique))
    for (i in 1:length(labels.unique)) {
        tt <- c("Cluster ", i)
        tt <- paste(tt, collapse = "")
        text[i] <- tt
    }
    lty <- rep(1, length(text))
    legend(x = min(x), y = max(x), legend = text, col = colours.templates.last, 
        lty = lty, cex = 0.6)
    sim.orig <- Result$similarity.orig
    sim.final <- Result$similarity.final
    etichette <- rep(0, 2)
    etichette[1] <- "orig. data"
    etichette[2] <- paste("k =", n.clust)
    dev.new()
    if (similarity.method == "d1.pearson" || similarity.method == 
        "d0.pearson" || similarity.method == "d1.pearson.mean" || 
        similarity.method == "d0.pearson.mean") {
        boxplot(sim.orig, sim.final, notch = FALSE, boxwex = 0.3, 
            col = c("grey", "orange"), ylim = c(min(sim.orig, 
                sim.final), 1))
    }
    if (similarity.method == "d0.L2" || similarity.method == 
        "d1.L2" || similarity.method == "d0.L2.centered" || similarity.method == 
        "d1.L2.centered") {
        boxplot(sim.orig, sim.final, notch = FALSE, boxwex = 0.3, 
            col = c("grey", "orange"), ylim = c(min(sim.orig, 
                sim.final), max(sim.orig, sim.final)))
    }
    title4 <- c("Registration: ", Result$warping.method)
    title4 <- paste(title4, collapse = "")
    title44 <- c("Boxplot Similarity Indexes")
    title4def <- c(title4, title44)
    title(main = title4def)
    axis(1, at = 1:2, labels = etichette, las = 0)
}
