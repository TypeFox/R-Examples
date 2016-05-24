aPvalResampOrd <- function(q, K, starResid)
## Evaluates the distribution function of a quadratic form in normal variables 
## at the value point 'q' following a standard permutation method. 'score' and 
## 'K' are respectively the vector and the matrix composing the quadratic form, 
## and 'starResid' is a matrix of permuted residuals. 
{
    inter1 <- starResid %*% K
    inter2 <- inter1 * starResid
    starQ <- apply(inter2, 1, sum)

    return(mean(starQ > q))
}
