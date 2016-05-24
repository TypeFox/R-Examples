pchi <-
function (q, df) 
{
    names.q <- names(q)
    arg.mat <- cbind.no.warn(q = as.vector(q), df = as.vector(df))
    na.index <- is.na.matrix(arg.mat)
    if (all(na.index)) 
        p <- rep(NA, nrow(arg.mat))
    else {
        p <- numeric(nrow(arg.mat))
        p[na.index] <- NA
        p.no.na <- p[!na.index]
        for (i in c("q", "df")) assign(i, arg.mat[!na.index, 
            i])
        if (any(df < .Machine$double.eps)) 
            stop("All non-missing values of 'df' must be positive")
        q.out <- q <= 0
        p.no.na[q.out] <- 0
        index <- !q.out
        if (any(index)) {
            p.no.na[index] <- pchisq(q[index]^2, df[index])
        }
        p[!na.index] <- p.no.na
    }
    if (!is.null(names.q)) 
        names(p) <- rep(names.q, length = length(p))
    else names(p) <- NULL
    p
}
