"equispaced" <-
function(y, m){
    if((n <- length(y)) != length(m))
        m <- rep(m, n)
    int <- lapply(1:n, function(i) seq(0, 1, length = m[i] + 2))
    a <- b <- numeric(n)
    for(i in 1:n){
        int. <- int[[i]]
        index <- findInterval(y[i] / m[i], int.)
        a[i] <- int.[index]
        b[i] <- int.[index + 1]
        if(is.na(b[i])){
            a[i] <- int.[m[i] + 1]
            b[i] <- 1
        }
    }
    cbind(a, b)
}

