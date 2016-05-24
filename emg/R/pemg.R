# This analytical form was contributed by Mark Kozdoba
# N(mu,sigma^2) + Exp(lambda)
pemg <- function (q, mu = 0, sigma = 1, lambda = 1, lower.tail = TRUE, log.p = FALSE)
{
    #renormalize
    u <- (q - mu) * lambda
    sigma1 <- lambda * sigma

    #compute for the case where mu=0, sigma = sigma1, lambda = 1
    p <- pnorm(u, 0, sigma1) - 
           exp(-u + (sigma1 * sigma1)/2 + 
               pnorm(u, sigma1 * sigma1, sigma1, log.p=TRUE) )
    
    p <- if(lower.tail) {p} else {1-p}
    
    if(log.p) {log(p)} else {p}
}