fitasy <-
function (y, B, b, p, c0) 
{
    a <- B %*% b
    ccr <- 0 * p
    for (j in 1:length(p)) {
        w <- ifelse(y >= c0[j] * a, p[j], 1 - p[j])
        ccr[j] <- sum(w * a * y)/sum(w * a * a)
    }
    return(ccr)
}
