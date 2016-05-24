
require(grDevices) # for gray

library(gridGraphics)

par1 <- function() {
    par("ylog") # FALSE
    plot(1 : 12, log = "y")
    par("ylog") # TRUE
}

par2 <- function() {
    plot(1:2, xaxs = "i") # 'inner axis' w/o extra space
    par(c("usr", "xaxp"))
}

par3 <- function() {
    ( nr.prof <-
     c(prof.pilots = 16, lawyers = 11, farmers = 10, salesmen = 9,
       physicians = 9, mechanics = 6, policemen = 6, managers = 6,
       engineers = 5, teachers = 4, housewives = 3, students = 3,
       armed.forces = 1))
    par(las = 3)
    barplot(rbind(nr.prof)) # R 0.63.2: shows alignment problem
    par(las = 0)  # reset to default
}

par4 <- function() {
    ## 'fg' use:
    plot(1:12, type = "b", main = "'fg' : axes, ticks and box in gray",
         fg = gray(0.7), bty = "7" , sub = R.version.string)
}

## Line types
showLty <- function(ltys, xoff = 0, ...) {
   stopifnot((n <- length(ltys)) >= 1)
   op <- par(mar = rep(.5,4)); on.exit(par(op))
   plot(0:1, 0:1, type = "n", axes = FALSE, ann = FALSE)
   y <- (n:1)/(n+1)
   clty <- as.character(ltys)
   mytext <- function(x, y, txt)
      text(x, y, txt, adj = c(0, -.3), cex = 0.8, ...)
   abline(h = y, lty = ltys, ...); mytext(xoff, y, clty)
   y <- y - 1/(3*(n+1))
   abline(h = y, lty = ltys, lwd = 2, ...)
   mytext(1/8+xoff, y, paste(clty," lwd = 2"))
}

par5 <- function() {
    showLty(c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash"))
    par(new = TRUE)  # the same:
    showLty(c("solid", "44", "13", "1343", "73", "2262"), xoff = .2, col = 2)
}

par6 <- function() {
    showLty(c("11", "22", "33", "44",   "12", "13", "14",   "21", "31"))
}

plotdiff(expression(par1()), "par-1")
plotdiff(expression(par2()), "par-2")
plotdiff(expression(par3()), "par-3")
plotdiff(expression(par4()), "par-4")
plotdiff(expression(par5()), "par-5")
plotdiff(expression(par6()), "par-6")

plotdiffResult()

