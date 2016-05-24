"describe2"<-function (x)
{
    lab <- c("N", "mean", "Std.Dev.","skew","kurt", "min", "Q1", "median",
        "Q3", "max", "missing values")
    if (missing(x)) {
        return(lab)
    }
    temp <- rep(0, length(lab))
    xt <- x[!is.na(x)]
    n <- length(xt)
    if (!is.numeric(xt) || all(is.na(x))) {
        return(c(n, NA, NA, rep(NA, 7), length(x) - length(xt)))
    }
    if (n < 4) {
        xt <- sort(xt)
        if (n == 1) {
            return(c(n, xt[1], NA, rep(xt[1], 7), length(x) -
                length(xt)))
        }
        if (n == 2) {
            return(c(n, mean(xt), sqrt(var(xt)),NA,NA, c(xt[1], xt[1],
                mean(xt), xt[2], xt[2]), length(x) - length(xt)))
        }
        if (n == 3) {
            return(c(n, mean(xt), sqrt(var(xt)),NA,NA, c(xt[1], xt[1],
                xt[2], xt[3], xt[3]), length(x) - length(xt)))
        }
    }
    else {
        return(c(length(xt), mean(xt), sqrt(var(xt)),skew(xt),kurt(xt), min(xt),
            quantile(xt, c(0.25, 0.5, 0.75)), max(xt), length(x) -
                length(xt)))
    }
}
