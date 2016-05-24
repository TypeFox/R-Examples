fmtcheck <-
function (x) 
{
    f1 <- (is.na(x)) * 1
    f2 <- (nchar(x) == 0) * 1
    g <- rep(0, length(x))
    for (i in 0:9) {
        ff <- (str_detect(x, as.character(i))) * 1
        g <- g + ff
    }
    f4 <- (g == 0) * 1
    f5 <- f1 + f2 + f4
}
