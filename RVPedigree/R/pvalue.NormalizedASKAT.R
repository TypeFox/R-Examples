##' Compute p-value and score for the normalized ASKAT method
##'
##' @inheritParams pvalue.ASKAT
##' @return A list with score and p-value for the Normalized ASKAT
##'     test on the given data. The list contains the elements
##'     \code{score} and \code{p.value}
##' @author Karim Oualkacha
##' @keywords internal
pvalue.NormalizedASKAT <- function(RH0, y, X, Phi, W, G){
    n        <- length(y)
    gam.hat  <- solve(t(X) %*% X) %*% t(X) %*% y
    residus  <- as.vector(y-X %*% gam.hat)
    eps      <- qnorm(rank(residus) / (n + 1))

    SS.ASKAT <- pvalue.ASKAT(RH0=RH0, y=eps, X=X, Phi=Phi, W=W, G=G)
    p.value  <- SS.ASKAT$p.value
    score    <- SS.ASKAT$score

    return(list(score=score, p.value=p.value))
}
