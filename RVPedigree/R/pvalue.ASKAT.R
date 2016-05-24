##' Compute p-value and score for the ASKAT method
##'
##' @param RH0 a vector of length 2 which the results (output) of the
##' \code{\link{Estim.H0.ASKAT}} function (i.e. variance components
##' estimates under the null model)
##' @inheritParams RVPedigree
##' @param G matrix of genotypes
##' @return A list with score and p-value for the ASKAT test on the
##'     given data. The list contains the elements \code{score} and
##'     \code{p.value}.
##' @author Karim Oualkacha
##' @importFrom CompQuadForm davies
##' @keywords internal
pvalue.ASKAT <- function(RH0, y, X, Phi, W, G)
{
    n      <- length(y)
    s.e    <- RH0[1]
    s.g    <- RH0[2]
    V      <- s.g * Phi + s.e * diag(rep(1, n))
    V.inv  <- solve(V)
    V.eig  <- eigen(V)
    V.sqrt <- V.eig$vectors %*% diag(sqrt(V.eig$values)) %*% solve(V.eig$vectors)
    P      <- V.inv -
        V.inv %*% X %*% solve(t(X) %*% V.inv %*% X) %*% t(X) %*% V.inv
    K      <- G %*% W %*% t(G)
    M      <- P %*% K %*% P
    Q      <- V.sqrt %*% M %*% V.sqrt

    Q.value <- eigen(Q)$values
    Q.value <- Q.value[Q.value > 0]

    score <- as.vector((y) %*% M %*% y)
    p.value <- davies(score, Q.value)$Qq
    return(list(score=score, p.value=p.value))
}
