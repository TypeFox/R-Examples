draw.lines <-
function (col = "red", alpha = 0.5, ...) 
{
  
    xy <- locator(2)
    lines(xy$x, xy$y, col = col.alpha(col[1], alpha), ...)

    N <- length(col)
    i <- 1
    while(is.list(c1 <- locator(1))) {
        i <- i + 1
        xy$x <- c(xy$x, c1$x)
        xy$y <- c(xy$y, c1$y)
        NL <- length(xy$x)
        lines(xy$x[(NL-1):NL], xy$y[(NL-1):NL], col = col.alpha(col[ifelse(i%%N == 0, N, i%%N)], alpha), ...)
        
    }
    invisible(data.frame(xy))
}

