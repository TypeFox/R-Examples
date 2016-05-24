##' Compute score and p-value for the VC-C1 method
##'
##' @param Prep a list of results obtained from the function
##' \code{\link{Preparation.VCC}}
##' @inheritParams RVPedigree
##' @param G genotype matrix
##' @return a list with the score and p-value of the association test
##'     in the given region. The list contains the elements
##'     \code{score} and \code{p.value}.
##' @author Karim Oualkacha
##' @author M'Hamed Lajmi Lakhal-Chaieb
##' @seealso \code{\link{pvalue.VCC2}}, \code{\link{pvalue.VCC3}}
##' @importFrom Matrix nearPD
##' @importFrom CompQuadForm davies
##' @keywords internal
pvalue.VCC1 <- function(Prep, W, G)
{
    L     <- Prep$L
    B     <- -Prep$B
    score <- t(L) %*% (G %*% W %*% t(G)) %*% L

    V <- nearPD(t(G) %*% B %*% G)$mat # ne pas avoir erreur (avoid
                                      # non-positive definte)
    V <- Matrix::as.matrix(V)
    V.eig     <- eigen(V)
    V.valeurs <- V.eig$values
    V.sqrt    <- V.eig$vectors %*% diag(sqrt(V.valeurs)) %*% solve(V.eig$vectors)

    psi.valeurs <- eigen(V.sqrt %*% W %*% V.sqrt)$values
    psi.valeurs <- psi.valeurs[psi.valeurs > 0]

    p.value <- davies(score, psi.valeurs)$Qq
    return(list(score=score, p.value=p.value))
}
