
require(stats)

library(gridGraphics)

set.seed(753)
bxp.data <- split(rt(100, 4), gl(5, 20))
bx.p <- boxplot(bxp.data, plot=FALSE)

bxp1 <- function() {
    boxplot(bxp.data)
}

bxp2 <- function() {
    par(mfrow =  c(2, 2))
    bxp(bx.p, xaxt = "n")
    bxp(bx.p, notch = TRUE, axes = FALSE, pch = 4, boxfill = 1:5)
    bxp(bx.p, notch = TRUE, boxfill = "lightblue", frame = FALSE,
        outl = FALSE, main = "bxp(*, frame= FALSE, outl= FALSE)")
    bxp(bx.p, notch = TRUE, boxfill = "lightblue", border = 2:6,
        ylim = c(-4,4), pch = 22, bg = "green", log = "x",
        main = "... log = 'x', ylim = *")
}
bxp3 <- function() {
    par(mfrow = c(1, 2))
    ## single group -- no label
    boxplot (weight ~ group, data = PlantGrowth, subset = group == "ctrl")
    ## with label
    bx <- boxplot(weight ~ group, data = PlantGrowth,
                  subset = group == "ctrl", plot = FALSE)
    bxp(bx, show.names=TRUE)
}

set.seed(1)
z <- split(rnorm(1000), rpois(1000, 2.2))

bxp4 <- function() {
    ## examples for new (S+ like) features
    boxplot(z, whisklty = 3, main = "boxplot(z, whisklty = 3)")
}

bxp5 <- function() {
    ## Colour support similar to plot.default:
    par(mfrow = 1:2, bg = "light gray", fg = "midnight blue")
    boxplot(z,   col.axis = "skyblue3",
            main = "boxplot(*, col.axis=..,main=..)")
    plot(z[[1]], col.axis = "skyblue3",
         main =    "plot(*, col.axis=..,main=..)")
    mtext("par(bg=\"light gray\", fg=\"midnight blue\")",
          outer = TRUE, line = -1.2)
}

bxp6 <- function() {
    ## Mimic S-Plus:
    splus <- list(boxwex = 0.4, staplewex = 1, outwex = 1, boxfill = "grey40",
                  medlwd = 3, medcol = "white", whisklty = 3, outlty = 1,
                  outpch = NA)
    boxplot(z, pars = splus)
    ## Recycled and "sweeping" parameters
}

bxp7 <- function() {
    par(mfrow = c(1,2))
    boxplot(z, border = 1:5, lty = 3, medlty = 1, medlwd = 2.5)
    boxplot(z, boxfill = 1:3, pch = 1:5, lwd = 1.5, medcol = "white")
}

bxp8 <- function() {
    ## too many possibilities
    boxplot(z, boxfill = "light gray", outpch = 21:25, outlty = 2,
            bg = "pink", lwd = 2,
            medcol = "dark blue", medcex = 2, medpch = 20)
}

plotdiff(expression(bxp1()), "bxp-1")
plotdiff(expression(bxp2()), "bxp-2")
plotdiff(expression(bxp3()), "bxp-3")
plotdiff(expression(bxp4()), "bxp-4")
plotdiff(expression(bxp5()), "bxp-5", width=8)
plotdiff(expression(bxp6()), "bxp-6")
plotdiff(expression(bxp7()), "bxp-7")
plotdiff(expression(bxp8()), "bxp-8")

plotdiffResult()
