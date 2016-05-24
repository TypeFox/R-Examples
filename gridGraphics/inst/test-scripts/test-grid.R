
library(gridGraphics)

grid1 <- function() {
    plot(1:3)
    grid(NA, 5, lwd = 2) # grid only in y-direction
}

grid2 <- function() {
    ## maybe change the desired number of tick marks:  par(lab = c(mx, my, 7))
    par(mfcol = 1:2)
    with(iris,
         {
             plot(Sepal.Length, Sepal.Width, col = as.integer(Species),
                  xlim = c(4, 8), ylim = c(2, 4.5), panel.first = grid(),
                  main = "with(iris,  plot(...., panel.first = grid(), ..) )")
             plot(Sepal.Length, Sepal.Width, col = as.integer(Species),
                  panel.first = grid(3, lty = 1, lwd = 2),
                  main = "... panel.first = grid(3, lty = 1, lwd = 2), ..")
         }
         )
}

plotdiff(expression(grid1()), "grid-1")
plotdiff(expression(grid2()), "grid-2", width=9)

plotdiffResult()
