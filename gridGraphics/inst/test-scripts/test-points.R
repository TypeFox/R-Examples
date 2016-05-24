
library(gridGraphics)

## ------------ test code for various pch specifications -------------
# Try this in various font families (including Hershey)
# and locales.  Use sign = -1 asserts we want Latin-1.
# Standard cases in a MBCS locale will not plot the top half.
TestChars <- function(sign = 1, font = 1, ...)
{
   MB <- l10n_info()$MBCS
   r <- if(font == 5) { sign <- 1; c(32:126, 160:254)
       } else if(MB) 32:126 else 32:255
   if (sign == -1) r <- c(32:126, 160:255)
   par(pty = "s")
   plot(c(-1,16), c(-1,16), type = "n", xlab = "", ylab = "",
        xaxs = "i", yaxs = "i",
        main = sprintf("sign = %d, font = %d", sign, font))
   grid(17, 17, lty = 1) ; mtext(paste("MBCS:", MB))
   for(i in r) try(points(i%%16, i%/%16, pch = sign*i, font = font,...))
}

points1 <- function() {
    TestChars()
}

points2 <- function() {
    try(TestChars(sign = -1))
}

points3 <- function() {
    TestChars(font = 5)  # Euro might be at 160 (0+10*16).
                         # Mac OS has apple at 240 (0+15*16).
}

points4 <- function() {
    try(TestChars(-1, font = 2))  # bold
}

plotdiff(expression(points1()), "points-1")
plotdiff(expression(points2()), "points-2")
plotdiff(expression(points3()), "points-3")
plotdiff(expression(points4()), "points-4")

plotdiffResult()
