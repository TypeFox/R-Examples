#' Update for total variation parameter in \code{equi_mcmc}.
#'
#' Samples from the square root of an inverse-gamma.
#'
#' This function provides a Gibbs update for the total variation parameter from
#' the MCMC implemented in \code{equi_mcmc}. This corresponds to the square root
#' of an inverse-gamma distributed random variable whose parameters depend on
#' the data and the component covariance matrices. Roughly, this is the update
#' for the standard deviation, not the variance.
#'
#' @param X An array. The tensor data.
#' @param phi_inv A list of the current values of inverse of the
#'   lower-triangular Cholesky square root of the the component covariance
#'   matrices. This is equivalent to the transpose of the upper-triangular
#'   Cholesky square root of the inverse component covariance matrices.
#'
#'   \code{phi_inv[[i]]} is a lower triangluar matrix where
#'   \code{solve(phi_inv[[i]]) \%*\% t(solve(phi_inv[[i]]))} is the current
#'   estimate of the \eqn{i}th component covariance matrix.
#' @return A numeric. The update for the total variation parameter in the MCMC
#'   implemented in \code{equi_bayes}.
#'
#' @seealso \code{\link{equi_mcmc}} for a Gibbs sampler where this function is
#'   used.
#'
#' @keywords equivariance
#'
#' @references Gerard, D., & Hoff, P. (2015).
#'   \href{http://www.sciencedirect.com/science/article/pii/S0047259X15000330}{Equivariant
#'    minimax dominators of the MLE in the array normal model}. \emph{Journal of
#'   Multivariate Analysis}, 137, 32-49.
#'
#' @author David Gerard.
sample_sig <- function(X, phi_inv) {
    p <- dim(X)
    alpha <- prod(p) / 2
    beta <- fnorm(atrans(X, phi_inv)) ^ 2 / 2
    new_sig <- 1 / sqrt(stats::rgamma(1, alpha, beta))
    return(new_sig)
}

#' Sample from the mirror-Wishart distribution.
#'
#' Given scale matrix \code{Phi} and degrees of freedom \code{nu},
#' \code{rmirror_wishart} will sample from the mirror-Wishart distribution.
#'
#' \eqn{S} is mirror-Wishart(\eqn{\nu,\Phi}) if \deqn{S = UV'VU',} where
#' \eqn{VV'} is the lower triangular Cholesky decomposition of a
#' Wishart(\eqn{\nu,I})-distributed random matrix and \eqn{UU'} is the upper
#' triangular Cholesky decomposition of \eqn{\Phi}. That is, \eqn{V} is lower
#' triangular and \eqn{U} is upper triangular. For details on its applications,
#' see
#' \href{http://www.sciencedirect.com/science/article/pii/S0047259X15000330}{
#' Gerard and Hoff (2015)}.
#'
#' @param nu An integer. The degrees of freedom in the mirror-Wishart.
#' @param Phi A matrix. The scale matrix of the mirror-Wishart.
#'
#' @return A matrix drawn from the mirror-Wishart distribution with \code{nu}
#'   degrees of freedom and scale matrix \code{Phi}.
#'
#' @export
#'
#' @references Gerard, D., & Hoff, P. (2015).
#'   \href{http://www.sciencedirect.com/science/article/pii/S0047259X15000330}{Equivariant
#'    minimax dominators of the MLE in the array normal model}. \emph{Journal of
#'   Multivariate Analysis}, 137, 32-49.
#'
#' @author David Gerard.
#'
#' @keywords equivariance simulation
#'
#' @seealso \code{\link{sample_right_wishart}}
rmirror_wishart <- function(nu, Phi) {
    n <- nrow(Phi)
    if (nu < n) {
        warning("nu must be >= n")
    }
    V <- matrix(0, nrow = n, ncol = n)
    V[lower.tri(V)] <- stats::rnorm(n * (n - 1) / 2)
    diag(V) <- sqrt(stats::rchisq(n, df = (nu - 1:n + 1)))

    #inefficient hack for upper triangular cholesky.
    U <- backsolve(chol(solve(Phi)), diag(n))

    return(U %*% t(V) %*% V %*% t(U))
}


#' Gibbs update of \code{Phi_inv}.
#'
#' Samples an upper triangular Cholesky square root of a
#' mirror-Wishart distributed random variable.
#'
#' Let \eqn{X} be mirror-Wishart(\eqn{\nu}, \eqn{V^-1}). Then This code
#' returns an upper triangular \eqn{C} where \eqn{X = CC'}. This
#' function is used primarily during the Gibbs updates of the inverse
#' of the lower triangular Cholesky square root of the component
#' covariance matrices in \code{equi_mcmc}.
#'
#' @param nu A numeric. The degrees of freedom in the mirror-Wishart.
#' @param V A matrix. The inverse of the scale matrix in the
#'     mirror-Wishart.
#'
#' @return \code{C} An upper triangular matrix such that \code{C \%*\% t(C)} is
#'   a sample from the mirror-Wishart(\code{nu}, \code{V ^ -1}) distribution.
#'
#' @seealso \code{\link{equi_mcmc}}, \code{\link{rmirror_wishart}}.
#'
#' @author David Gerard.
#'
#' @keywords equivariance simulation
#'
#' @references Gerard, D., & Hoff, P. (2015).
#'   \href{http://www.sciencedirect.com/science/article/pii/S0047259X15000330}{Equivariant
#'   minimax dominators of the MLE in the array normal model}. \emph{Journal of
#'   Multivariate Analysis}, 137, 32-49.
#'
sample_right_wishart <- function(nu, V) {
    m <- nrow(V)
    df <- (nu + nu - m + 1) - (nu - m + 1):nu
    if (m > 1) {
        Tmat <- diag(sqrt(stats::rchisq(c(rep(1, m)), df)))
        Tmat[lower.tri(Tmat)] <- stats::rnorm(m * (m + 1) / 2 - m)
    } else {
        Tmat <- sqrt(stats::rchisq(1, df))
    }
    U <- chol(V)
    C <- backsolve(U, diag(m)) %*% t(Tmat)
    return(C)
}

#' Gibbs sampler using an invariant prior.
#'
#' \code{equi_mcmc} obtains posterior draws that are useful in optimal
#' equivariant estimation under the array normal model.
#'
#' \code{equi_mcmc} obtains posterior samples of the component
#' covariance matrices from the array normal model. This is with
#' respect to using the right Haar measure over a product group of
#' lower triangular matrices as the prior.
#'
#' This returns only the upper triangular Cholesky square root of the
#' inverses of the component covariance matrices. Equivalently, these
#' are the inverses of the lower triangular Cholesky square roots of
#' the component covariance matrices. This is because sampling the
#' inverse is faster computationally and the Bayes rules (based on
#' multiway Stein's loss) only depend on the inverse.
#'
#' @param X A tensor.
#' @param itermax The number of iterations in the Gibb's sampler.
#' @param start_identity Should we start the component covariance
#'     matrices at the identity (TRUE) or the sample covariance
#'     matrices (FALSE)?
#' @param print_iter Should we print the iteration number at each
#'     iteration?
#' @param mode_rep The mode that contains samples. I.e., the mode
#'     whose component covariance matrix is the identity. If NULL then
#'     no modes are assumed to have identity covariance.
#'
#' @return \code{Phi_inv} List of posterior draws of the inverse of
#'     the cholesky square roots of each component covariance
#'     matrix. \code{Phi_inv[[i]][,,j]} provides the \eqn{j}th sample
#'     of the \eqn{i}th component.
#'
#'   \code{sigma} Vector of posterior samples of the overall scale
#'   paramater.
#'
#' @seealso \code{\link{sample_right_wishart}} and
#'     \code{\link{sample_sig}} for the Gibbs
#'     updates. \code{\link{convert_cov}} and
#'     \code{\link{get_equi_bayes}} for getting posterior summaries
#'     based on the output of
#'     \code{equi_mcmc}. \code{\link{multiway_takemura}} for an
#'     improvement on this procedure.
#'
#' @export
#'
#' @keywords equivariance
#'
#' @author David Gerard.
#'
#' @references Gerard, D., & Hoff, P. (2015).
#'   \href{http://www.sciencedirect.com/science/article/pii/S0047259X15000330}{Equivariant
#'   minimax dominators of the MLE in the array normal model}. \emph{Journal of
#'   Multivariate Analysis}, 137, 32-49.
#'
#' @examples
#' #Generate data whose true covariance is just the identity.
#' p <- c(2,2,2)
#' X <- array(stats::rnorm(prod(p)),dim = p)
#' #Then run the Gibbs sampler.
#' mcmc_out <- equi_mcmc(X)
#' plot(mcmc_out$sigma, type = 'l', lwd = 2, ylab = expression(sigma),
#'      xlab = 'Iteration', main = 'Trace Plot')
#' abline(h = 1,col = 2,lty = 2)
equi_mcmc <- function(X, itermax = 1000, start_identity = FALSE, print_iter = FALSE, mode_rep = NULL) {
    p <- dim(X)
    n <- length(p)
    phi_inv <- list()
    sigma <- rep(NA, length = itermax)
    sigma_solo <- rep(NA, length = itermax)

    ## Initiate parameter values
    sigma[1] <- 1
    sigma_solo[1] <- 1
    if (start_identity) {
        ## then inverse starts at identity
        start_ident_list <- start_ident(p)
        for (mode_index in 1:n) {
            phi_inv[[mode_index]] <- array(NA, dim = c(p[mode_index], p[mode_index], itermax))
            phi_inv[[mode_index]][, , 1] <- start_ident_list[[mode_index]]
        }
    } else {
        ## then inverse starts at the inverse of the residual sums of squares
        start_resid_list <- start_resids(X)
        for (mode_index in 1:n) {
            phi_inv[[mode_index]] <- array(NA, dim = c(p[mode_index], p[mode_index], itermax))
            phi_inv[[mode_index]][, , 1] <- t(backsolve(t(start_resid_list[[2]][[mode_index]]), diag(p[mode_index])))
        }
    }

    if (!is.null(mode_rep)) {
        phi_inv[[mode_rep]][, , 1] <- diag(p[mode_rep])
    } else {
        mode_rep <- n + 1
    }


    for (iter_index in 2:itermax) {
        phi_inv_current <- list()
        for (mode_index in 1:n) {
            phi_inv_current[[mode_index]] <- phi_inv[[mode_index]][, , iter_index - 1]
        }
        sigma_solo[iter_index] <- sample_sig(X, phi_inv_current)
        sigma[iter_index] <- sigma_solo[iter_index]

        for (mode_index in (1:n)[-mode_rep]) {
            nu <- prod(p[-mode_index])
            phi_mult <- phi_inv_current
            phi_mult[[mode_index]] <- diag(p[mode_index])
            half_V <- mat(atrans(X, phi_mult), mode_index)
            V <- tcrossprod(half_V, half_V)
            L_temp <- t(sample_right_wishart(nu, V))
            sigma[iter_index] <- det(L_temp) ^ (-1 / p[mode_index])
            phi_inv[[mode_index]][, , iter_index] <- L_temp * sigma[iter_index]
            phi_inv_current[[mode_index]] <- phi_inv[[mode_index]][, , iter_index]
        }
        if (!identical(mode_rep, n + 1)) {
            phi_inv[[mode_rep]][, , iter_index] <- diag(p[mode_rep])
            phi_inv_current[[mode_rep]] <- phi_inv[[mode_rep]][, , iter_index]
        }


        if (print_iter) {
            cat("Iteration =", iter_index, "\n")
        }
    }
    return(list(Phi_inv = phi_inv, sigma = sigma))
}
