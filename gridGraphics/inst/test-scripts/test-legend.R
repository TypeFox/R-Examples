
library(gridGraphics)

require(stats)

legend1 <- function() {
    ## Run the example in '?matplot' or the following:
    leg.txt <- c("Setosa     Petals", "Setosa     Sepals",
                 "Versicolor Petals", "Versicolor Sepals")
    y.leg <- c(4.5, 3, 2.1, 1.4, .7)
    cexv  <- c(1.2, 1, 4/5, 2/3, 1/2)
    matplot(c(1, 8), c(0, 4.5), type = "n", xlab = "Length", ylab = "Width",
            main = "Petal and Sepal Dimensions in Iris Blossoms")
    for (i in seq(cexv)) {
        text  (1, y.leg[i] - 0.1, paste("cex=", formatC(cexv[i])),
               cex = 0.8, adj = 0)
        legend(3, y.leg[i], leg.txt, pch = "sSvV", col = c(1, 3), cex = cexv[i])
    }
}

legend2 <- function() {
    ## 'merge = TRUE' for merging lines & points:
    x <- seq(-pi, pi, len = 65)
    plot(x, sin(x), type = "l", ylim = c(-1.2, 1.8), col = 3, lty = 2)
    points(x, cos(x), pch = 3, col = 4)
    lines(x, tan(x), type = "b", lty = 1, pch = 4, col = 6)
    title("legend(..., lty = c(2, -1, 1), pch = c(NA, 3, 4), merge = TRUE)",
          cex.main = 1.1)
    legend(-1, 1.9, c("sin", "cos", "tan"), col = c(3, 4, 6),
           text.col = "green4", lty = c(2, -1, 1), pch = c(NA, 3, 4),
           merge = TRUE, bg = "gray90")
}

legend3 <- function() {
    ## right-justifying a set of labels: thanks to Uwe Ligges
    x <- 1:5; y1 <- 1/x; y2 <- 2/x
    plot(rep(x, 2), c(y1, y2), type = "n", xlab = "x", ylab = "y")
    lines(x, y1); lines(x, y2, lty = 2)
    temp <- legend("topright", legend = c(" ", " "),
                   text.width = strwidth("1,000,000"),
                   lty = 1:2, xjust = 1, yjust = 1,
                   title = "Line Types")
    text(temp$rect$left + temp$rect$w, temp$text$y,
         c("1,000", "1,000,000"), pos = 2)
}

legend4 <- function() {
    ##--- log scaled Examples ------------------------------
    leg.txt <- c("a one", "a two")

    par(mfrow = c(2, 2))
    for(ll in c("","x","y","xy")) {
        plot(2:10, log = ll, main = paste0("log = '", ll, "'"))
        abline(1, 1)
        lines(2:3, 3:4, col = 2)
        points(2, 2, col = 3)
        rect(2, 3, 3, 2, col = 4)
        text(c(3,3), 2:3, c("rect(2,3,3,2, col=4)",
                            "text(c(3,3),2:3,\"c(rect(...)\")"),
             adj = c(0, 0.3))
        legend(list(x = 2,y = 8), legend = leg.txt, col = 2:3, pch = 1:2,
               lty = 1, merge = TRUE)   #, trace = TRUE)
    }
}

legend5 <- function() {
    ##-- Math expressions:  ------------------------------
    x <- seq(-pi, pi, len = 65)
    plot(x, sin(x), type = "l", col = 2, xlab = expression(phi),
         ylab = expression(f(phi)))
    abline(h = -1:1, v = pi/2*(-6:6), col = "gray90")
    lines(x, cos(x), col = 3, lty = 2)
    ex.cs1 <- expression(plain(sin) * phi,  paste("cos", phi))  # 2 ways
    utils::str(legend(-3, .9, ex.cs1, lty = 1:2, plot = FALSE,
                      adj = c(0, 0.6)))  # adj y !
    legend(-3, 0.9, ex.cs1, lty = 1:2, col = 2:3,  adj = c(0, 0.6))
}

legend6 <- function() {
    x <- rexp(100, rate = .5)
    hist(x, main = "Mean and Median of a Skewed Distribution")
    abline(v = mean(x),   col = 2, lty = 2, lwd = 2)
    abline(v = median(x), col = 3, lty = 3, lwd = 2)
    ex12 <- expression(bar(x) == sum(over(x[i], n), i == 1, n),
                       hat(x) == median(x[i], i == 1, n))
    utils::str(legend(4.1, 30, ex12, col = 2:3, lty = 2:3, lwd = 2))
}

legend7 <- function() {
    ## 'Filled' boxes -- for more, see example(plot.factor)
    op <- par(bg = "white") # to get an opaque box for the legend
    plot(cut(weight, 3) ~ group, data = PlantGrowth, col = NULL,
         density = 16*(1:3))
}

legend8 <- function() {
    ## Using 'ncol' :
    x <- 0:64/64
    matplot(x, outer(x, 1:7, function(x, k) sin(k * pi * x)),
            type = "o", col = 1:7, ylim = c(-1, 1.5), pch = "*")
    op <- par(bg = "antiquewhite1")
    legend(0, 1.5, paste("sin(", 1:7, "pi * x)"), col = 1:7, lty = 1:7,
           pch = "*", ncol = 4, cex = 0.8)
    legend(.8,1.2, paste("sin(", 1:7, "pi * x)"), col = 1:7, lty = 1:7,
           pch = "*", cex = 0.8)
    legend(0, -.1, paste("sin(", 1:4, "pi * x)"), col = 1:4, lty = 1:4,
           ncol = 2, cex = 0.8)
    legend(0, -.4, paste("sin(", 5:7, "pi * x)"), col = 4:6,  pch = 24,
           ncol = 2, cex = 1.5, lwd = 2, pt.bg = "pink", pt.cex = 1:3)
}

legend9 <- function() {
    ## point covering line :
    x <- 0:64/64
    y <- sin(3*pi*x)
    plot(x, y, type = "l", col = "blue",
         main = "points with bg & legend(*, pt.bg)")
    points(x, y, pch = 21, bg = "white")
    legend(.4,1, "sin(c x)", pch = 21, pt.bg = "white", lty = 1, col = "blue")
}

legend10 <- function() {
    ## legends with titles at different locations
    x <- 0:64/64
    y <- sin(3*pi*x)
    plot(x, y, type = "n")
    legend("bottomright", "(x,y)", pch = 1, title = "bottomright")
    legend("bottom", "(x,y)", pch = 1, title = "bottom")
    legend("bottomleft", "(x,y)", pch = 1, title = "bottomleft")
    legend("left", "(x,y)", pch = 1, title = "left")
    legend("topleft", "(x,y)", pch = 1, title = "topleft, inset = .05",
           inset = .05)
    legend("top", "(x,y)", pch = 1, title = "top")
    legend("topright", "(x,y)", pch = 1, title = "topright, inset = .02",
           inset = .02)
    legend("right", "(x,y)", pch = 1, title = "right")
    legend("center", "(x,y)", pch = 1, title = "center")
}

legend11 <- function() {
    # using text.font (and text.col):
    par(mfrow = c(2, 2), mar = rep(2.1, 4))
    c6 <- terrain.colors(10)[1:6]
    for(i in 1:4) {
        plot(1, type = "n", axes = FALSE, ann = FALSE);
        title(paste("text.font =",i))
        legend("top", legend = LETTERS[1:6], col = c6,
               ncol = 2, cex = 2, lwd = 3, text.font = i, text.col = c6)
    }
}

plotdiff(expression(legend1()), "legend-1")
plotdiff(expression(legend2()), "legend-2", density=72)
plotdiff(expression(legend3()), "legend-3")
plotdiff(expression(legend4()), "legend-4")
plotdiff(expression(legend5()), "legend-5")
plotdiff(expression(legend6()), "legend-6")
plotdiff(expression(legend7()), "legend-7")
plotdiff(expression(legend8()), "legend-8")
plotdiff(expression(legend9()), "legend-9")
plotdiff(expression(legend10()), "legend-10")
plotdiff(expression(legend11()), "legend-11")

plotdiffResult()
