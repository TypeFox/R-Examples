dnhyper <- function(x, m, n, k, log = FALSE){
    if (k > m) {
        stop(paste("k cannot be larger than m!",
            "\n"))
    }
    p <- exp(lchoose(x - 1, k - 1) + lchoose(n - x, m - k) - lchoose(n, m))
    if (log)
        p <- log(p)
    p[is.nan(p)] <- 0
    p <- pmin(pmax(p, 0),1)
    p
}

pnhyper <- function(q, m, n, k, lower.tail = TRUE, log.p = FALSE){
    if (k > m) {
        stop(paste("k cannot be larger than m!",
            "\n"))
    }
    q <- ceiling(q)
    temp <- rep(0, length(q))
    for (i in 1:length(q)) {
        if (q[i] < k) {
            temp[i] <- 0
        }
        else if (q[i] >= (n - m + k)) {
            temp[i] <- 1
        }
        else temp[i] <- sum(dnhyper(k:q[i], m = m, n = n, k = k))
    }
    if (lower.tail == FALSE)
        temp <- 1 - temp
    if (log.p)
        temp <- log(temp)
    temp <- pmin(pmax(temp, 0),1)
    temp
}

qnhyper <- function(p, m, n, k, lower.tail = TRUE, log.p = FALSE){
    if (k > m) {
        stop(paste("k cannot be larger than m!",
            "\n"))
    }
    if (log.p) p <- exp(p)
    all.p <- NULL
    temp <- k:(n - m + k)
    if(lower.tail){
        temp.out <- rbind(c(k,-Inf), cbind(temp, pnhyper(temp, m, n, k)), c(n-m+k,Inf))
    } else temp.out <- rbind(c(k,Inf),cbind(temp, pnhyper(temp, m, n, k, lower.tail = FALSE)),c(n-m+k,-Inf))
    for (i in 1:length(p)) {
    if(lower.tail){
        all.p <- c(all.p, min(temp.out[which(temp.out[,2]>=p[i]),1]))
    } else all.p <- c(all.p, min(temp.out[which(temp.out[,2]<p[i]),1]))
       }
    all.p
}

rnhyper <- function(nn, m, n, k){
    if (k > m) {
        stop(paste("k cannot be larger than m!",
            "\n"))
    }
    qnhyper(runif(nn), m, n, k)
}