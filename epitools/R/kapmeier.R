"kapmeier" <-
function (time, status) 
{
    stime <- sort(time)
    status <- status[order(time)]
    nj <- length(time):1
    nj <- nj[!duplicated(stime)]
    dj <- tapply(status, stime, sum)
    tj <- unique(stime)
    sj <- (nj - dj)/nj
    cumsj <- cumprod(sj)
    cumrj <- 1 - cumsj
    results <- cbind(time = tj, n.risk = nj, n.events = dj, condsurv = sj, 
        survival = cumsj, risk = cumrj)
    dimnames(results)[1] <- list(NULL)
    results[dj != 0, ]
}

