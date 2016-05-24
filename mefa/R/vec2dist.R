vec2dist <-
function(x, size, labels = NULL, diag = FALSE, upper = FALSE, call = FALSE, method = NULL)
{
    class(x) <- "dist"
    attr(x, "Size") <- size
    attr(x, "Labels") <- labels
    attr(x, "Diag") <- diag
    attr(x, "Upper") <- upper
    attr(x, "call") <- if (call) match.call() else NULL
    attr(x, "method") <- method
    x
}
