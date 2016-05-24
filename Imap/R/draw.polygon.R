draw.polygon <-
function (col = "blue", alpha = 0.5, lty = 1, ...) 
{
    col <- col.alpha(col, alpha)

    xy <- locator(2)
    lines(xy$x, xy$y, lty = lty)

    while(is.list(c1 <- locator(1))) {
        xy$x <- c(xy$x, c1$x)
        xy$y <- c(xy$y, c1$y)
        lines(xy$x, xy$y, lty = lty)
    }
    xy <- data.frame(xy)
    xy <- rbind(xy, xy[1, ])
    polygon(xy$x, xy$y, lty = lty, col = col, ...)
    
    invisible(xy)
}

