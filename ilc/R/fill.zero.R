fill.zero <- function (x, method = "constant") 
{
    tt <- 1:length(x)
    zeros <- abs(x) < 1e-09
    xx <- x[!zeros]
    tt <- tt[!zeros]
    x <- approx(tt, xx, 1:length(x), method = method, f = 0.5, 
        rule = 2)
    return(x$y)
}
