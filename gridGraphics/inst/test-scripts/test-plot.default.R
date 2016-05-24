
library(gridGraphics)

Speed <- cars$speed
Distance <- cars$dist

plot.default1 <- function() {
    plot(Speed, Distance, panel.first = grid(8, 8),
         pch = 0, cex = 1.2, col = "blue")
}

plot.default2 <- function() {
    plot(Speed, Distance,
         panel.first = lines(stats::lowess(Speed, Distance), lty = "dashed"),
         pch = 0, cex = 1.2, col = "blue")
}

plot.default3 <- function() {
    ## Show the different plot types
    x <- 0:12
    y <- sin(pi/5 * x)
    par(mfrow = c(3,3), mar = .1+ c(2,2,3,1))
    for (tp in c("p","l","b",  "c","o","h",  "s","S","n")) {
        plot(y ~ x, type = tp, main = paste0("plot(*, type = \"", tp, "\")"))
        if(tp == "S") {
            lines(x, y, type = "s", col = "red", lty = 2)
            mtext("lines(*, type = \"s\", ...)", col = "red", cex = 0.8)
        }
    }
}

plot.default4 <- function() {
    ##--- Log-Log Plot  with  custom axes
    lx <- seq(1, 5, length = 41)
    yl <- expression(e^{-frac(1,2) * {log[10](x)}^2})
    y <- exp(-.5*lx^2)
    par(mfrow = c(2,1), mar = c(5, 5, 1, 1))
    plot(10^lx, y, log = "xy", type = "l", col = "purple",
         main = "Log-Log plot", ylab = yl, xlab = "x")
    plot(10^lx, y, log = "xy", type = "o", pch = ".", col = "forestgreen",
         main = "Log-Log plot with custom axes", ylab = yl, xlab = "x",
         axes = FALSE, frame.plot = TRUE)
    my.at <- 10^(1:5)
    axis(1, at = my.at, labels = formatC(my.at, format = "fg"))
    at.y <- 10^(-5:-1)
    axis(2, at = at.y, labels = formatC(at.y, format = "fg"), col.axis = "red")
}

plotdiff(expression(plot.default1()), "plot.default-1")
plotdiff(expression(plot.default2()), "plot.default-2")
plotdiff(expression(plot.default3()), "plot.default-3")
plotdiff(expression(plot.default4()), "plot.default-4", height=12)

plotdiffResult()
