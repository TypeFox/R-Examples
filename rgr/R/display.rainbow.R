display.rainbow <-
function () 
{
    frame()
    palette(rainbow(36))
    x <- 0:7
    y <- 0:7
    plot(x, y, xaxt = "n", yaxt = "n", xlab = "", ylab = "", 
        type = "n", main = "The rainbow(36) palette")
    for (i in 1:6) for (j in 1:6) {
        ii <- (i - 1) * 6 + j
        text(i, j, paste(ii), col = ii, adj = 0.5, cex = 2)
    }
    palette("default")
    invisible()
}
