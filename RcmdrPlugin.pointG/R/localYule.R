localYule<-
function (X) 
{
    sumrow <- apply(X, 1, sum)
    sumcol <- apply(X, 2, sum)
    n <- sum(X)

    result <- X
    for (i in (1:nrow(X))) {
    for (j in (1:ncol(X))) {
        a <- X[i, j]
        b <- sumrow[i] - a
        c <- sumcol[j] - a
        d <- n - a - b - c
        Q <- (a * d - b * c)/(a * d + b * c)
        result[i,j] <- Q
    }
}
    result
}

