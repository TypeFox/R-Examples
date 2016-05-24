
library(gridGraphics)

require(stats); require(grDevices)

set.seed(1)
x <- 1:10
y <- sort(10*runif(10))
z <- runif(10)
z3 <- cbind(z, 2*runif(10), runif(10))

symbols1 <- function() {
    symbols(x, y, thermometers = cbind(.5, 1, z), inches = .5, fg = 1:10)
}

symbols2 <- function() {
    symbols(x, y, thermometers = z3, inches = FALSE)
    text(x, y, apply(format(round(z3, digits = 2)), 1, paste, collapse = ","),
         adj = c(-.2,0), cex = .75, col = "purple", xpd = NA)
}

## Note that  example(trees)  shows more sensible plots!
N <- nrow(trees)

symbols3 <- function() {
    with(trees, {
         ## Girth is diameter in inches
         symbols(Height, Volume, circles = Girth/24, inches = FALSE,
                 main = "Trees' Girth") # xlab and ylab automatically
     })
}

symbols4 <- function() {
    ## Colours too:
    palette(rainbow(N, end = 0.9))
    with(trees, {
        symbols(Height, Volume, circles = Girth/16, inches = FALSE, bg = 1:N,
                fg = "gray30",
                main = "symbols(*, circles = Girth/16, bg = 1:N)")
    })
}

# Some of my own tests to cover the range of symbols
symbols5 <- function() {
    symbols(mtcars$disp, mtcars$mpg,
            rect=as.matrix(abs(scale(mtcars[, c(3, 1)]))))
}

symbols6 <- function() {
    symbols(mtcars$disp, mtcars$mpg,
            stars=as.matrix(abs(scale(mtcars[, c(4:7, 10:11)]))))
}

symbols7 <- function() {
    x <- y <- w <- h <- lw <- uw <- 1:5
    m <- 1:5/6
    symbols(x, y, boxplots=cbind(w, h, lw, uw, m))
}

plotdiff(expression(symbols1()), "symbols-1")
plotdiff(expression(symbols2()), "symbols-2")
plotdiff(expression(symbols3()), "symbols-3")
plotdiff(expression(symbols4()), "symbols-4", density=72, antialias=FALSE)
plotdiff(expression(symbols3()), "symbols-5")
plotdiff(expression(symbols3()), "symbols-6")
plotdiff(expression(symbols3()), "symbols-7")

plotdiffResult()
