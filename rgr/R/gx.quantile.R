gx.quantile <-
function (xx, f, display = TRUE) 
{
    temp.x <- remove.na(xx)
    x <- sort(temp.x$x[1:temp.x$n])
    pr <- (rank(x) - 0.5)/temp.x$n
    q <- signif(approx(pr, x, f)$y, 4)
    if (display) 
        cat("  Fractile for", deparse(substitute(xx)), "=", f, 
            "\n  Quantile =", q, "\n")
    invisible(list(q = q))
}
