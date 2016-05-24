aPvalResampMom <- function(q, K, starResid)
## Evaluates the distribution function of a quadratic form in normal variables 
## at the value point 'q' following a permutation method based on matching 
## moments. 'score' and 'K' are respectively the vector and the matrix composing
## the quadratic form, and 'starResid' is a matrix of permuted residuals.
{
    inter1 <- starResid %*% K
    inter2 <- inter1 * starResid
    starQ <- apply(inter2, 1, sum)
    
    B <- length(starQ)
    muQ <- mean(starQ)
    sig2Q <- (B - 1) / B * var(starQ)
    mu4 <- sum((starQ - muQ)^4) / B
    gamma <- mu4 / sig2Q^2 - 3
    dfd <- 12 / gamma
    return(pchisq((q - muQ) * sqrt(2 * dfd) / sqrt(sig2Q) + dfd, df=dfd, 
                  lower.tail=FALSE))
}
