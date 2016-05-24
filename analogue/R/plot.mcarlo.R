plot.mcarlo <-
    function (x, which = c(1:2), alpha = 0.05,
              caption = c("Distribution of dissimilarities",
              expression(paste("Simulated probability Pr (Dissim < ", alpha,
                  ")"))),
              col.poly = "lightgrey", border.poly = "lightgrey",
              ask = prod(par("mfcol")) < length(which) && dev.interactive(),
              ...)
{
    if (!inherits(x, "mcarlo"))
        stop("\"x\" must be or inherit from class \"mcarlo\".")
    if (!is.numeric(which) || any(which < 1) || any(which > 2))
        stop("'which' must be in 1:2")
    show <- rep(FALSE, 2)
    show[which] <- TRUE
    one.fig <- prod(par("mfcol")) == 1
    if (ask) {
        op <- par(ask = TRUE)
        on.exit(par(op))
    }
    xlabel <- paste("Dissimilarity (", attr(x, "method"), ")",
                    sep = "")
    evalDists <- switch(attr(x, "method"),
                        SQchord = seq(from = 0, to = 2, by = 0.01),
                        chord = seq(from = 0, to = sqrt(2), by = 0.01),
                        seq(from = 0, to = max(x, na.rm = TRUE), by = 0.01))
    if (show[1]) {
        range.dists <- max(evalDists)
        h <- hist(x, plot = FALSE)
        d <- density(x, adjust = 2, from = 0, to = range.dists)
        y.lim <- range(0, h$density, d$y)
        hist(x, freq = FALSE, ylim = y.lim, xlim = c(0, range.dists),
             main = "", xlab = "", ylab = "", ann = FALSE, col = col.poly,
             border = border.poly)
        box()
        title(main = caption[1], ylab = "Density", xlab = xlabel)
        lines(d)
    }
    if (show[2]) {
        num.dists <- length(x)
        cumfreq <- sapply(evalDists,
                          function(x, y) length(y[y <= x]), x,
                          USE.NAMES = FALSE) / num.dists
        cummu <- data.frame(distances = evalDists, cumfreq = cumfreq)
        suppressWarnings(critical <- approx(cummu$cumfreq, cummu$distances,
                                            alpha))
        tlen <- length(critical$x)
        sigx <- cummu$distances <= (max(critical$y, na.rm = TRUE))
        sigx <- cummu[sigx, ]
        plot(cumfreq ~ distances, data = cummu, type = "n", ann = FALSE)
        title(main = caption[2], ylab = "Probability", xlab = xlabel)
        polygon(c(0, sigx[, 1], sigx[nrow(sigx), 1]),
                c(0, sigx[, 2], 0), col = col.poly, border = border.poly)
        lines(cumfreq ~ distances, data = cummu)
    }
    if(show[2])
        invisible(critical)
    else
        invisible()
}
