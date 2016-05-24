rblnorm <-
function(n, prob, meanlog, sdlog)
{
    if(max(length(prob), length(meanlog), length(sdlog)) > 1)
        stop("parameters must be of length 1")
    p <- runif(n)
    q <- rep(0, length(p))
    cases <- p > (1-prob)
    q[cases] <- rlnorm(sum(cases), meanlog=meanlog, sdlog=sdlog)
    q
}

