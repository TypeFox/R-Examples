rbgamma <-
function(n, prob, scale, shape)
{
    if(max(length(prob), length(scale), length(shape)) > 1)
        stop("parameters must be of length 1")
    p <- runif(n)
    q <- rep(0, length(p))
    cases <- p > (1-prob)
    q[cases] <- rgamma(sum(cases), scale=scale, shape=shape)
    q
}

