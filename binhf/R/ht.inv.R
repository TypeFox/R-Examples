"ht.inv" <-
function (data) 
{
    if (is.list(data)) {
        nsc <- length(data$s)
        data <- c(data$d, data$s[nsc])
    }
    a <- 2
    n <- length(data)
    nhalf <- n/2
    J <- logb(n, 2)
    res <- sm1 <- NULL
    for (i in 1:J) {
        res <- c(data[1:nhalf],res)
        data <- data[(nhalf + 1):n]
        n <- n/2
        nhalf <- nhalf/2
    }
    res <- c(data[1], res)
    nhalf <- 1
    n <- 2
    sm <- rep(0, nhalf)
    det <- sm
    for (i in 1:J) {
        sm[1:nhalf] <- res[1:nhalf]
        sm1 <- c(sm[1:nhalf], sm1)
        det[1:nhalf] <- res[(nhalf + 1):n]
        res[2 * (1:nhalf) - 1] <- a/2 * (sm[1:nhalf] + det[1:nhalf])
        res[2 * (1:nhalf)] <- a/2 * (sm[1:nhalf] - det[1:nhalf])
        n <- 2 * n
        nhalf <- 2 * nhalf
    }
    return(list(res = res, sm1 = sm1))
}

