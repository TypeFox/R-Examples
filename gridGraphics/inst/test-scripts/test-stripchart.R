
library(gridGraphics)

set.seed(1)
x <- stats::rnorm(50)
xr <- round(x, 1)

stripchart1 <- function() {
    stripchart(x) ; m <- mean(par("usr")[1:2])
    text(m, 1.04, "stripchart(x, \"overplot\")")
    stripchart(xr, method = "stack", add = TRUE, at = 1.2)
    text(m, 1.35, "stripchart(round(x,1), \"stack\")")
    stripchart(xr, method = "jitter", add = TRUE, at = 0.7)
    text(m, 0.85, "stripchart(round(x,1), \"jitter\")")
}

stripchart2 <- function() {
    stripchart(decrease ~ treatment,
               main = "stripchart(OrchardSprays)",
               vertical = TRUE, log = "y", data = OrchardSprays)
}

stripchart3 <- function() {
    stripchart(decrease ~ treatment, at = c(1:8)^2,
               main = "stripchart(OrchardSprays)",
               vertical = TRUE, log = "y", data = OrchardSprays)
}

plotdiff(expression(stripchart1()), "stripchart-1")
plotdiff(expression(stripchart2()), "stripchart-2")
plotdiff(expression(stripchart3()), "stripchart-3")

plotdiffResult()

