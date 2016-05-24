##' Prepare for the VC-C methods
##'
##' @param n sample size (number of individuals)
##' @param RH0 a list of the results obtained from the function
##' \code{\link{Estim.H0.VCC}}
##' @param Phi kinship matrix
##' @return XXX
##' @author Karim Oualkacha
##' @author M'Hamed Lajmi Lakhal-Chaieb
##' @importFrom ks kdde
##' @importFrom ks hns
##' @keywords internal
Preparation.VCC <- function(n, RH0, Phi)
{
    h <- RH0$h
    gam <- RH0$gam
    residus <- RH0$residuals

    F <- rank(residus) / (n + 1)
    S <- qnorm(F)

    h0 <- hns(residus, deriv.order=0)
    f0 <- kdde(residus, h=h0, deriv.order=0, eval.points=residus)$estimate

    h1 <- hns(residus, deriv.order=1)
    f1 <- kdde(residus, h=h1, deriv.order=1, eval.points=residus)$estimate

    h2 <- hns(residus, deriv.order=2)
    f2 <- kdde(residus, h=h2, deriv.order=2, eval.points=residus)$estimate

    F.deriv1 <- 1 / dnorm(qnorm(F))
    F.deriv2 <- qnorm(F) * F.deriv1^2

    s1 <- F.deriv1 * f0
    s2 <- F.deriv2 * (f0^2) + F.deriv1 * f1

    tilde.s1 <- diag(s1)
    tilde.s2 <- diag(s2)

    K <- log(f0)
    k1 <- f1 / f0
    k2 <- (f2 / f0) - ((f1 / f0)^2)

    tilde.k2 <- diag(k2)

    I.N <- diag(rep(1, n))
    Gamma <- h * Phi + (1 - h) * I.N
    V <- I.N - solve(Gamma)

    tilde.VS <- diag(as.vector(V %*% S))

    L <- k1 + tilde.s1 %*% (V %*% S)
    B <- tilde.s1 %*% V %*% tilde.s1 +
        tilde.s2 %*% tilde.VS + tilde.k2

    P <- B + L %*% t(L)

    return(list(P=P, L=L, B=B))
}
