qqLine <-
function (x, y = NULL, ...) 
{
    x.quartiles <- quantile(x, c(0.25, 0.75))
    if (is.null(y)) {
        norm.quartiles <- qnorm(c(0.25, 0.75))
        b <- (x.quartiles[2] - x.quartiles[1])/(norm.quartiles[2] - 
            norm.quartiles[1])
        a <- x.quartiles[1] - norm.quartiles[1] * b
    }
    else {
        y.quartiles <- quantile(y, c(0.25, 0.75))
        b <- (y.quartiles[2] - y.quartiles[1])/(x.quartiles[2] - 
            x.quartiles[1])
        a <- y.quartiles[1] - x.quartiles[1] * b
    }
    abline(a, b, ...)
    ans <- as.vector(c(a, b))
    names(ans) <- c("intercept", "slope")
    invisible(ans)
}
