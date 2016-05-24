##' Estimate the model parameters under the null model
##'
##' This function estimates the model parameters under the null model
##' when there is no region genotypes effect, for the VCC methods
##' @inheritParams Estim.H0.ASKAT
##' @return a list of:
##' \itemize{
##' \item \code{h.hat}: estimate of the (narrow sense) heritability
##' \item \code{gam.hat}: a vector of estimates of fixed effects of the covariates
##' (if there are any).
##' \item \code{residuals}: residuals from adjusting the null model
##' with only covariates in the model.
##' }
##' @author Karim Oualkacha
##' @author M'Hamed Lajmi Lakhal-Chaieb
##' @export
Estim.H0.VCC <- function(y, X, S, U)
{
    y <- check_pheno(y)
    check_covariates(X, y)

    gam.hat <- solve(t(X) %*% X) %*% t(X) %*% y
    residus <- as.vector(y - X %*% gam.hat)
    n       <- length(y)
    eps     <- qnorm(rank(residus) / (n + 1))
    h.hat   <- optimize(Neg.LogLikelihood.VC.C1,
                        interval=c(0, 1),
                        eps=eps,
                        U=U,
                        S=S)$min

    return(list(h=h.hat, gam=gam.hat, residuals=residus))
}
