#' Get list of identity matrices.
#'
#' Will provide a list of identity matrices for the specified modes.
#'
#' Given a vector of dimensions \code{p} and a vector indicating which
#' modes will get an identity matrix \code{modes}, this function will
#' return a list \code{start_vals} where \code{start_vals[[i]]} is the
#' identity matrix of dimensions \code{p[i]} if \code{i} is in
#' \code{modes} and will be \code{NULL} otherwise.
#'
#' This is primarily used when getting starting values in \code{equi_mcmc}.
#'
#' @param p A vector of integers. This is the dimension of the array and the
#'   length of the list to be created.
#' @param modes A vector of integers. These are the indices in the list to be
#'   given an identity matrix.
#'
#' @return \code{start_vals} A list of identity matrices and NULL values.
#'
#' @export
#'
#' @author David Gerard.
#'
#' @seealso \code{\link{equi_mcmc}}.
start_ident <- function(p, modes = NULL) {
    n <- length(p)
    if (is.null(modes)) {
        modes <- 1:n
    }
    start_vals <- list()
    for (index in modes) {
        start_vals[[index]] <- diag(p[index])
    }
    return(start_vals)
}

#' Sample covariance matrices for each mode.
#'
#' Scaled Cholesky square roots of the sample covariance matrix and
#' its inverse.
#'
#' This function will take the sample covariance matrix of the
#' \eqn{i}th matricization of an input array \eqn{Y} and will return
#' (1) its lower-triangular Cholesky square root scaled down to have
#' determinant 1 and (2) the inverse of its lower-triangular Cholesky
#' square root scaled down to have determinant 1. This function is
#' primarily used to obtain starting values for the Gibbs sampler
#' implemented in \code{equi_mcmc}.
#'
#' @param Y An array of numeric data.
#' @param mode_rep A vector of integers. The modes specified by
#'     \code{mode_rep} will be given an identity matrix instead of a
#'     sample-based matrix.
#' @return \code{Sig} A list where \code{Sig[[i]]} is the
#'     lower-triangular Cholesky square root of the sample covariance
#'     matrix of the \eqn{i}th mode, scaled down to have determinant
#'     1.
#'
#'   \code{Sig_inv} A list where \code{Sig_inv[[i]]} is the inverse of the
#'   lower-triangular Cholesky square root of the sample covariance matrix of
#'   the \eqn{i}th mode, scaled down to have determinant 1.
#'
#'   If \code{mode_rep} is not \code{NULL}, then the list elements in \code{Sig}
#'   and \code{Sig_inv} specified in \code{mode_rep} will be the identity matrix
#'   instead of sample-based matrices.
#'
#' @export
#'
#' @author David Gerard.
#'
#' @seealso \code{\link{equi_mcmc}}.
start_resids <- function(Y, mode_rep = NULL) {
    p <- dim(Y)
    K <- length(p)
    Sig <- list()
    Sig_inv <- list()
    if (is.null(mode_rep)) {
        mode_rep <- K + 1
    }
    for (k in (1:K)[-mode_rep]) {
        Sig[[k]] <- mat(Y, k) %*% t(mat(Y, k))
        Sig[[k]] <- t(chol(Sig[[k]] / det(Sig[[k]]) ^ (1 / p[k])))
        Sig_inv[[k]] <- backsolve(Sig[[k]], diag(p[k]), upper.tri = FALSE)
    }
    if (max(mode_rep) <= K) {
        for (k in mode_rep) {
            if (k <= K) {
                Sig[[k]] <- diag(p[k])
                Sig_inv[[k]] <- diag(p[k])
            }
        }
    }
    return(list(Sig = Sig, Sig_inv = Sig_inv))
}

#' Demeans array data.
#'
#' Rotates an array into two parts, one of which has mean zero.
#'
#' If one mode contains samples (or repetitions), then this function
#' will rotate the array into two parts, a mean part and a covariance
#' part. The 'covariance part' has mean zero and the rest of the
#' methods in this package apply. The 'mean part' is simply the sample
#' mean. If the data are array normal, then the 'covariance part' will
#' also be array normal with the \emph{exact} same covariance
#' structure as the original tensor, except that there are one fewer
#' samples.
#'
#' @param X An array, one of whose modes is assumed to be samples from
#'     the array normal model.
#' @param mode_reps The mode(s) that contain(s) the samples, or
#'     repetitions, from the array normal model.
#'
#' @author David Gerard.
#'
#' @references Gerard, D., & Hoff, P. (2015).
#'   \href{http://www.sciencedirect.com/science/article/pii/S0047259X15000330}{Equivariant
#'   minimax dominators of the MLE in the array normal model}. \emph{Journal of
#'   Multivariate Analysis}, 137, 32-49.
#'
#' @return \code{Y} An array that has the same dimensions as \code{X}
#'     except that the mode \code{mode_reps} has dimension one
#'     smaller. This array is mean 0 array normal with the same
#'     covariance structure as \code{X}.
#'
#'   \code{X_bar} The sample mean of \code{X}. Under the array normal
#'   model, \code{X} and \code{Y} are statistically independent.
#'
#' @export
demean_tensor <- function(X, mode_reps) {
    p <- dim(X)
    eigen_vecs <- (eigen(diag(p[mode_reps]) -
                           matrix(1 / p[mode_reps],nrow = p[mode_reps], ncol = p[mode_reps]))$vectors)
    H <- eigen_vecs[, -p[mode_reps]]
    Y <- amprod(X, t(H), mode_reps)
    # X_bar <- amprod(X,t(eigen_vecs[,p[mode_reps]]),mode_reps) / sqrt(p[mode_reps]) ## equivalent calculation to MLE
    X_bar <- apply(X, (1:length(p))[-mode_reps], mean)

    return(list(Y = Y, X_bar = X_bar))
}
