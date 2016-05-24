myprettyW <-
function (x, extradigits = 0)
{
    if (is.na(x))
        return(NA)
    if (x < 1) {
        n <- 3
    }
    else if (x > 1000) {
        n <- ceiling(log(x, 10))
    }
    else {
        n <- 4
    }
    n <- n + extradigits
    signif(x, n)
}
