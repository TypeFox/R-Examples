
library(gridGraphics)

boxplot1 <- function() {
    ## boxplot on a formula:
    boxplot(count ~ spray, data = InsectSprays, col = "lightgray")
    # *add* notches (somewhat funny here):
    boxplot(count ~ spray, data = InsectSprays,
            notch = TRUE, add = TRUE, col = "blue")
}

boxplot2 <- function() {
    boxplot(decrease ~ treatment, data = OrchardSprays,
            log = "y", col = "bisque")
}

boxplot3 <- function() {
    rb <- boxplot(decrease ~ treatment, data = OrchardSprays, col = "bisque")
    title("Comparing boxplot()s and non-robust mean +/- SD")
    mn.t <- tapply(OrchardSprays$decrease, OrchardSprays$treatment, mean)
    sd.t <- tapply(OrchardSprays$decrease, OrchardSprays$treatment, sd)
    xi <- 0.3 + seq(rb$n)
    points(xi, mn.t, col = "orange", pch = 18)
    arrows(xi, mn.t - sd.t, xi, mn.t + sd.t,
           code = 3, col = "pink", angle = 75, length = .1)
}

mat <- cbind(Uni05 = (1:100)/21, Norm = rnorm(100),
             `5T` = rt(100, df = 5), Gam2 = rgamma(100, shape = 2))

boxplot4 <- function() {
    ## boxplot on a matrix:
    boxplot(as.data.frame(mat),
            main = "boxplot(as.data.frame(mat), main = ...)")
}

boxplot5 <- function() {
    par(las = 1) # all axis labels horizontal
    boxplot(as.data.frame(mat), main = "boxplot(*, horizontal = TRUE)",
            horizontal = TRUE)
}

boxplot6 <- function() {
    ## Using 'at = ' and adding boxplots -- example idea by Roger Bivand :
    boxplot(len ~ dose, data = ToothGrowth,
            boxwex = 0.25, at = 1:3 - 0.2,
            subset = supp == "VC", col = "yellow",
            main = "Guinea Pigs' Tooth Growth",
            xlab = "Vitamin C dose mg",
            ylab = "tooth length",
            xlim = c(0.5, 3.5), ylim = c(0, 35), yaxs = "i")
    boxplot(len ~ dose, data = ToothGrowth, add = TRUE,
            boxwex = 0.25, at = 1:3 + 0.2,
            subset = supp == "OJ", col = "orange")
    legend(2, 9, c("Ascorbic acid", "Orange juice"),
           fill = c("yellow", "orange"))
}

plotdiff(expression(boxplot1()), "boxplot-1")
plotdiff(expression(boxplot2()), "boxplot-2")
# Was getting weird single-pixel difference in (anti-aliased) dashed line
plotdiff(expression(boxplot3()), "boxplot-3", width=9, height=9)
plotdiff(expression(boxplot4()), "boxplot-4")
plotdiff(expression(boxplot5()), "boxplot-5")
plotdiff(expression(boxplot6()), "boxplot-6")

plotdiffResult()
