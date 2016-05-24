#' Estimation of the variance components under the null model using the ASKAT method
#'
#' @param y Vector of phenotype values
#' @param X A matrix of covariates, including intercept.
#' @param S Matrix obtained from spectral decomposition of the
#' relationship matrix: \eqn{\Phi = U S U^T}.
#' @param U Matrix obtained from spectral decomposition of the
#' relationship matrix: \eqn{\Phi = U S U^T}.
#' @return A vector with the following values:
#' \itemize{
#' \item \code{s.e}: variance component due to the error term (residuals)
#' \item \code{s.g}: variance component due to the polygenic residual
#' (polygenic background)
#' \item \code{beta}: the estimate of the effects of the covariate (if
#' there are some in the data)
#' }
#'
#' @author Karim Oualkacha
#' @author M'Hamed Lajmi Lakhal-Chaieb
#' @export
Estim.H0.ASKAT <- function(y, X, S, U)
{
    y <- check_pheno(y)
    check_covariates(X, y)

    Ut.x <- t(U) %*% X
    Ut.y <- t(U) %*% y
    delta <- optimize(Neg.LogLikelihood.ASKAT,
                      interval=c(1e-4, 100),
                      S=S,
                      Ut.y=Ut.y,
                      Ut.x=Ut.x,
                      n=length(y))$min

    W <- diag(1/(delta + S))
    beta <- solve(t(Ut.x) %*% W %*% Ut.x) %*% t(Ut.x) %*% W %*% Ut.y
    s.g <- mean((Ut.y - Ut.x %*% beta)^2 / (delta + S))
    s.e <- delta * s.g
    return(c(s.e, s.g, beta))
}
