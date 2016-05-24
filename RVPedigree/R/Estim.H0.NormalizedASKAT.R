##' Estimation of the variance components under the null model using
##' the normalized ASKAT method
##'
##' @inheritParams Estim.H0.ASKAT
##' @return A vector with the following values:
#' \itemize{
#' \item \code{s.e}: variance component due to the error term (residuals)
#' \item \code{s.g}: variance component due to the polygenic residual
#' (polygenic background)
#' \item \code{beta}: the estimate of the effects of the covariate (if
#' there are some in the data)
#' }
##' @author Karim Oualkacha
##' @author M'Hamed Lajmi Lakhal-Chaieb
##' @export
Estim.H0.NormalizedASKAT <- function(y, X, S, U)
{
    y <- check_pheno(y)
    check_covariates(X, y)

    gam.hat <- solve(t(X) %*% X) %*% t(X) %*% y
    residus <- as.vector(y - X %*% gam.hat)
    n       <- length(y)
    eps     <- qnorm(rank(residus) / (n + 1))

    Res.H0.SKAT <- Estim.H0.ASKAT(y = eps, X = X, S = S, U = U)

    return(Res.H0.SKAT)
}
