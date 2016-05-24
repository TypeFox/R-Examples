
library(gridGraphics)

## NASA space shuttle o-ring failures
fail <- factor(c(2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 2, 1, 2, 1, 1, 1,
                 1, 2, 1, 1, 1, 1, 1),
               levels = 1:2, labels = c("no", "yes"))
temperature <- c(53, 57, 58, 63, 66, 67, 67, 67, 68, 69, 70, 70,
                 70, 70, 72, 73, 75, 75, 76, 76, 78, 79, 81)

cdplot1 <- function() {
    ## CD plot
    cdplot(fail ~ temperature)
}
cdplot2 <- function() {
    cdplot(fail ~ temperature, bw = 2)
}

cdplot3 <- function() {
    cdplot(fail ~ temperature, bw = "SJ")
}

cdplot4 <- function() {
    ## compare with spinogram
    (spineplot(fail ~ temperature, breaks = 3))
}

cdplot5 <- function() {
    ## highlighting for failures
    cdplot(fail ~ temperature, ylevels = 2:1)
}

cdplot6 <- function() {
    ## scatter plot with conditional density
    cdens <- cdplot(fail ~ temperature, plot = FALSE)
    plot(I(as.numeric(fail) - 1) ~ jitter(temperature, factor = 2),
         xlab = "Temperature", ylab = "Conditional failure probability")
    lines(53:81, 1 - cdens[[1]](53:81), col = 2)
}

plotdiff(expression(cdplot1()), "cdplot-1")
plotdiff(expression(cdplot2()), "cdplot-2")
plotdiff(expression(cdplot3()), "cdplot-3")
plotdiff(expression(cdplot4()), "cdplot-4", width=9)
plotdiff(expression(cdplot5()), "cdplot-5")
plotdiff(expression(cdplot6()), "cdplot-6")

plotdiffResult()
