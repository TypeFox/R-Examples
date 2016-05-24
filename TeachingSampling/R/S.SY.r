S.SY<-function (N, a)
{
    r <- sample(a, 1)
    c <-  N - a * floor(N/a)
    if (r <= c)
        n <- floor((N/a)) + 1
    else n <- floor(N/a)
    sam <- matrix(0, n, 1)
    for (k in 0:n) {
        sam[k] <- r + (a * (k - 1))
    }
    sam
}

