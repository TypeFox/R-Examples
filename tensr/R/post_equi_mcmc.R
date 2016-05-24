#' Convert the output from \code{equi_mcmc} to component covariance matrices.
#'
#' This takes the output from \code{equi_mcmc}, which are the inverses of the
#' lower-triangular Cholesky square roots of the component covariance matrices,
#' and returns the component covariance matrices. These are the more useful
#' posterior draws to use in actual data analysis.
#'
#' The output from \code{equi_mcmc} is the inverse of the lower-triangular
#' Cholesky square root of each component covariance matrix. This output is
#' convenient for calculating the Bayes rule under multiway-Stein's loss (see
#' \code{\link{get_equi_bayes}}). Call one of these outputs from
#' \code{equi_mcmc} \eqn{\Psi}. Then this function calculates \eqn{\Sigma =
#' \Psi^-1\Psi^-T}, which are the posterior draws of the component covariance
#' matrices. These component covariance matrices are constrained to have
#' determinant one, hence there is a total variation parameter \eqn{\sigma^2}.
#'
#' @param equi_mcmc_obj The output from \code{equi_mcmc}, which contains a list.
#'   The first element is a list containing the posterior draws of the inverses
#'   of the lower-triangular Cholesky square roots of each component covariance
#'   matrix. The second list element is a total variation parameter, but the
#'   square root of the version used in calculating the overall covariance
#'   matrix.
#'
#' @return \code{cov_post} A list containing the posterior draws of each
#'   component covariance matrix.
#'
#'   \code{sig2_post} A vector containing the posterior draws of the total
#'   variation parameter.
#'
#' @seealso \code{\link{equi_mcmc}}.
#'
#' @export
#'
#' @keywords equivariance posterior
#'
#' @author David Gerard.
#'
#' @references  Gerard, D., & Hoff, P. (2015).
#'   \href{http://www.sciencedirect.com/science/article/pii/S0047259X15000330}{Equivariant
#'   minimax dominators of the MLE in the array normal model}. \emph{Journal of
#'   Multivariate Analysis}, 137, 32-49.
#'
#' @examples
#' #Generate data whose true covariance is just the identity.
#' p <- c(4,4,4)
#' X <- array(stats::rnorm(prod(p)),dim = p)
#' #Then run the Gibbs sampler.
#' mcmc_out <- equi_mcmc(X)
#' cov_out <- convert_cov(mcmc_out)
#'
#' # Some trace plots.
#' plot(cov_out[[2]], type = 'l', xlab = 'Iteration',
#'      ylab = expression(sigma ^ 2), main = 'Trace Plot')
#' abline(h = 1, col = 2, lty = 2)
#' legend('topleft', 'True Value', col = 2, lty = 2, bty = 'n')
#'
#' k <- sample(1:length(p), size = 1)
#' i <- sample(1:p[k], size = 1)
#' j <- sample(1:p[k], size = 1)
#' plot(cov_out[[1]][[k]][i, j, ], type = 'l', xlab = 'Iteration',
#'      main = 'Trace Plot',
#'      ylab = substitute(Sigma[k][group('[', list(i, j), ']')],
#'                        list(k = k, i = i, j = j)))
#' abline(h = 1 * (i == j), lty =  2, col = 2)
#' legend('topleft', 'True Value', col = 2, lty = 2, bty = 'n')
convert_cov <- function(equi_mcmc_obj) {
    ## This function returns posterior draws of each component covariance matrix.
    n <- length(equi_mcmc_obj[[1]])
    p <- rep(NA, length = n)
    ## convert to posterior draws of variance
    sig2_post <- equi_mcmc_obj[[2]] ^ 2


    chol_half_list <- list()
    for (mode_index in 1:n) {
        p[mode_index] <- dim(equi_mcmc_obj[[1]][[mode_index]])[1]
        chol_half_list[[mode_index]] <- array(apply(equi_mcmc_obj[[1]][[mode_index]], 3, backsolve,
                                                    upper.tri = FALSE, x = diag(p[mode_index])),
                                              dim = dim(equi_mcmc_obj[[1]][[mode_index]]))
    }
    cov_list <- list()
    for (mode_index in 1:n) {
        cov_list[[mode_index]] <- array(apply(chol_half_list[[mode_index]], 3, tcrossprod),
                                        dim = dim(equi_mcmc_obj[[1]][[mode_index]]))
    }
    return(list(cov_post = cov_list, sig2_post = sig2_post))
}

#' Get the Bayes rule under multiway Stein's loss.
#'
#' Given the output of \code{equi_mcmc}, this function will calculate the Bayes
#' rule under multiway Stein's loss.
#'
#' Multiway Stein's loss is a generalization of Stein's loss to more than two
#' dimensions. The Bayes rule under this loss is simply represented in terms of
#' the posterior moments of the component precision matrices. These moments can
#' be approximated by using the output of \code{equi_mcmc}. When using the
#' invariant prior that is used in \code{equi_mcmc}, the resulting Bayes rule is
#' the uniformly minimum risk equivariant estimator.
#'
#' More details on multiway Stein's loss and the Bayes rules under it can be
#' found in
#' \href{http://www.sciencedirect.com/science/article/pii/S0047259X15000330}{
#' Gerard and Hoff (2015)}.
#'
#' @param psi_inv A list of arrays where \code{psi_inv[[i]][[, , j]]} is the
#'   \eqn{j}th update of the \code{i}th component. These components are the
#'   inverses of the lower-triangular Cholesky square roots of the component
#'   covariance matrices. You can just use the \code{Phi_inv} output from
#'   \code{equi_mcmc}.
#' @param sigma A vector of posteior draws of the total variation parameter.
#'   This is just \code{sigma} from the output of \code{equi_mcmc}.
#' @param burnin A numeric between 0 and 1. What proportion of the posterior
#'   samples do you want to discard as burnin? The default is 0.25.
#' @return \code{Sig_hat} A list of the Bayes rules of the component covariance
#'   matrices under multiway Stein's loss.
#'
#'   \code{B} A list of the lower-triangular Cholesky square roots of the Bayes
#'   rules of the component covariance matrices under multiway Stein's loss. We
#'   have that \code{Sig_hat[[i]]} is equal to \code{B[[i]] \%*\% t(B[[i]])}.
#'
#'   \code{b} A numeric. This is the bayes rule of the total variation
#'   parameter. This is the 'standard deviation' version. That is, the \code{b ^
#'   2} would be used to calculate the overall covariance matrix.
#'
#' @seealso \code{\link{equi_mcmc}}.
#'
#' @export
#'
#' @keywords equivariance posterior
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
#' p <- c(4,4,4)
#' X <- array(stats::rnorm(prod(p)),dim = p)
#' #Then run the Gibbs sampler.
#' mcmc_out <- equi_mcmc(X)
#' bayes_rules <- get_equi_bayes(mcmc_out$Phi_inv, mcmc_out$sigma)
#' bayes_rules$Sig_hat[[1]]
get_equi_bayes <- function(psi_inv, sigma, burnin = NULL) {
    itermax <- length(sigma)
    if (is.null(burnin)) {
        burnin <- round(itermax * 0.25)
    }
    p <- sapply(psi_inv, dim)[1, ]
    n <- length(p)
    Sig_hat_list <- vector(mode = "list", length = n)
    ## initialize
    for (mode_index in 1:n) {
        Sig_hat_list[[mode_index]] <- t(psi_inv[[mode_index]][, , burnin]) %*%
          psi_inv[[mode_index]][, , burnin] / sigma[burnin] ^ 2
    }
    ## Now sum over all the iterations
    for (iter_index in (burnin + 1):itermax) {
        for (mode_index in 1:n) {
            Sig_hat_list[[mode_index]] <- Sig_hat_list[[mode_index]] + t(psi_inv[[mode_index]][, , iter_index]) %*%
              psi_inv[[mode_index]][, , iter_index] / sigma[iter_index] ^ 2
        }
    }
    for (mode_index in 1:n) {
        ## take inverse of average to get estimator
        Sig_hat_list[[mode_index]] <- solve(Sig_hat_list[[mode_index]] / (itermax - burnin + 1))
    }

    ## set the determinant to 1 and find scale estimate
    det_sig <- rep(NA, length = n)
    B_list <- vector(mode = "list", length = n)
    Sig_hat <- vector(mode = "list", length = n)
    b_sum <- 0
    for (mode_index in 1:n) {
        det_sig[mode_index] <- det(Sig_hat_list[[mode_index]])
        B_list[[mode_index]] <- t(chol(Sig_hat_list[[mode_index]]))
        B_list[[mode_index]] <- B_list[[mode_index]] / det_sig[mode_index] ^ (1 / (2 * p[mode_index]))
        Sig_hat[[mode_index]] <- B_list[[mode_index]] %*% t(B_list[[mode_index]])
        b_sum <- b_sum + det_sig[mode_index] ^ (1 / -p[mode_index])
    }
    b_hat <- n / b_sum
    return(list(B = B_list, b = sqrt(b_hat), Sig_hat = Sig_hat))
}


#' Generate a list of orthogonal matrices drawn from Haar distribution.
#'
#' Given a vector \code{p}, \code{random_ortho} will generate a list
#' \code{ortho_list} such that \code{ortho_list[[i]]} is a matrix with row and
#' column dimensions \code{p[[i]]} and is drawn from the uniform (Haar)
#' distribution over the space of orthogonal matrices.
#'
#' This function is primarily used by \code{\link{multiway_takemura}} in its
#' averaging over uniformly minimum risk equivariant estimators under rotations
#' of the data array.
#'
#' @param p A vector of dimensions for the matrices.
#'
#' @return \code{ortho_list} A list of orthogonal matrices whose dimensions are
#'   given in \code{p}.
#'
#' @export
#'
#' @keywords equivariance simulation
#'
#' @author David Gerard.
#'
#' @seealso \code{\link{multiway_takemura}}.
random_ortho <- function(p) {
    n <- length(p)
    ortho_list <- list()
    for (index in 1:n) {
        x_temp <- matrix(stats::rnorm(p[index] ^ 2), ncol = p[index], nrow = p[index])
        ortho_list[[index]] <- solve(mhalf(x_temp %*% t(x_temp))) %*% x_temp
    }
    return(ortho_list)
}

#' Calculate a truncated multiway Takemura estimator.
#'
#' This function will 'average' Bayes rules given random rotations of the data
#' array. This 'averaged' estimator has lower risk than the uniformly minimum
#' risk equivariant estimator under a product group of lower triangular
#' matrices. Truncated multiway Takemura's estimator is not equivariant with
#' respect to this product group of lower triangular matrices, but it is an
#' equivariant randomized estimator with respect to a product group of
#' orthogonal matrices.
#'
#' This function will (1) randomly rotate \code{X} along every mode, then (2) it
#' will calculate the uniformly minimum risk equivariant estimator using
#' \code{equi_mcmc}, then (3) it will 'average' these estimates.
#'
#' @param X An array. This is the data array.
#' @param ortho_max An integer. The number of 'averagings' to perform.
#' @param mcmc_itermax An integer. The number of iterations each MCMC should
#'   perform using \code{equi_mcmc}.
#' @param start_identity Should each MCMC start their covariance matrices at the
#'   identity (TRUE) or at the sample covariance matrices (FALSE)?
#' @param print_mcmc Should the output of the MCMC be printed to the screen
#'   (TRUE) or not (FALSE)?
#' @param mode_rep A vector of integers. Which mode(s) are considered iid
#'   observations? Default is none.
#'
#' @return \code{B} A list of the truncated multiway Takemura's estimators for
#'   each component covariance matrix. Not their Cholesky square roots.
#'
#'   \code{b} Truncated multiway Takemura's estimator for the total variation
#'   parameter. The 'variance' form, not the 'standard devation' form.
#'
#' @export
#'
#' @author David Gerard.
#'
#' @keywords equivariance
#'
#' @seealso \code{\link{equi_mcmc}}, \code{\link{random_ortho}}.
#'
#' @references Gerard, D., & Hoff, P. (2015).
#'   \href{http://www.sciencedirect.com/science/article/pii/S0047259X15000330}{Equivariant
#'   minimax dominators of the MLE in the array normal model}. \emph{Journal of
#'   Multivariate Analysis}, 137, 32-49.
#'
#' @examples
#' # Simulate data.
#' p <- c(5, 5, 5)
#' X <- array(stats::rnorm(prod(p)), dim = p)
#' multi_out <- multiway_takemura(X, mode_rep = 3)
#' multi_out$b
#' trim(multi_out$B[[1]])
#' trim(multi_out$B[[2]])
#' trim(multi_out$B[[3]])
multiway_takemura <- function(X, ortho_max = 2, mcmc_itermax = 1000, start_identity = FALSE,
                              print_mcmc = FALSE, mode_rep = NULL) {
    p <- dim(X)
    n <- length(p)
    for (ortho_index in 1:ortho_max) {
        ## transform orthogonally along each mode
        ortho_list <- random_ortho(p)
        X_ortho <- atrans(X, ortho_list)

        bayes_fit <- equi_mcmc(X_ortho, itermax = mcmc_itermax, start_identity = TRUE, mode_rep = mode_rep)
        bayes_rule_ortho <- get_equi_bayes(bayes_fit$Phi_inv, bayes_fit$sigma)

        ## orthogonallize back
        bayes_rule_back <- list()
        bayes_rule_back[[1]] <- list()
        bayes_rule_back[[2]] <- bayes_rule_ortho$b
        for (bayes_index in 1:n) {
            temp_bayes_rule <- t(ortho_list[[bayes_index]]) %*% bayes_rule_ortho$B[[bayes_index]]
            bayes_rule_back[[1]][[bayes_index]] <- temp_bayes_rule %*% t(temp_bayes_rule)
        }

        if (ortho_index == 1) {
            bayes_rule <- bayes_rule_back
        } else {
            ## now weighted average to get final bayes rule, but do the averaging on the trace scale,
            ## then scale back down to determinant scale.
            bayes_rule[[2]] <- bayes_rule_back[[2]] / ortho_index + bayes_rule[[2]] * (ortho_index - 1) / ortho_index
            names(bayes_rule) <- c("B", "b")
            for (bayes_index in 1:n) {
                back_temp <- bayes_rule_back[[1]][[bayes_index]]
                back_temp <- back_temp / tr(back_temp)
                bayes_prev_temp <- bayes_rule[[1]][[bayes_index]] / tr(bayes_rule[[1]][[bayes_index]])
                combo_temp <- back_temp / ortho_index + (ortho_index - 1) / ortho_index * bayes_prev_temp
                bayes_rule[[1]][[bayes_index]] <- combo_temp / det(combo_temp) ^ (1 / p[bayes_index])
            }
        }
        if (print_mcmc) {
            cat("Number of averagings = ", ortho_index, "\n")
        }
    }
    return(bayes_rule)
}
