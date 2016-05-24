dchi <-
function (x, df) 
{
    names.x <- names(x)
    arg.mat <- cbind.no.warn(x = as.vector(x), df = as.vector(df))
    na.index <- is.na.matrix(arg.mat)
    if (all(na.index)) 
        y <- rep(NA, nrow(arg.mat))
    else {
        y <- numeric(nrow(arg.mat))
        y[na.index] <- NA
        y.no.na <- y[!na.index]
        for (i in c("x", "df")) assign(i, arg.mat[!na.index, 
            i])
        if (any(df < .Machine$double.eps)) 
            stop("All non-missing values of 'df' must be positive")
        x.out <- x <= 0 | x == Inf
        y.no.na[x.out] <- 0
        index <- !x.out
        if (any(index)) {
            y.no.na[index] <- dchisq(x[index]^2, df[index]) * 
                2 * x[index]
        }
        y[!na.index] <- y.no.na
    }
    if (!is.null(names.x)) 
        names(y) <- rep(names.x, length = length(y))
    y
}
