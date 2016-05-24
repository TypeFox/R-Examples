
library(gridGraphics)

boxplot.matrix1 <- function() {
    set.seed(1)
    mat <- cbind(Uni05 = (1:100)/21, Norm = rnorm(100),
                 T5 = rt(100, df = 5), Gam2 = rgamma(100, shape = 2))
    boxplot(mat, main = "boxplot.matrix(...., main = ...)",
            notch = TRUE, col = 1:4)
}

plotdiff(expression(boxplot.matrix1()), "boxplot.matrix-1")

plotdiffResult()
