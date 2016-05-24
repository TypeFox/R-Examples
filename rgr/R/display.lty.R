display.lty <-
function () 
{
    frame()
    x <- 0:10
    y <- 0:10
    plot(x, y, xaxt = "n", yaxt = "n", xlab = "", ylab = "", 
        type = "n", main = "Available Line Types (lty) and Colours (colr)")
    for (i in 1:9) {
        ypos <- 10 - i
        text(1, ypos, paste(i), adj = 0.5, cex = 2)
        lines(x = c(3, 7), y = c(ypos, ypos), lty = i, col = 1)
        text(9, ypos, paste(i), col = i, adj = 0.5, cex = 2)
    }
    invisible()
}
