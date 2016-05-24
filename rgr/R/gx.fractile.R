gx.fractile <-
function (xx, q, display = TRUE) 
{
    temp.x <- remove.na(xx)
    x <- sort(temp.x$x[1:temp.x$n])
    n <- temp.x$n
    pr <- (rank(x) - 0.5)/n
    f <- signif((approx(x, pr, q)$y), 4)
    if (display) 
        cat("  Value of", deparse(substitute(xx)), "=", q, "\n  Fractile =", 
            f, "\n")
    invisible(list(f = f))
}
