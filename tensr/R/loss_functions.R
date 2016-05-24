
#' Calculate multiway Stein's loss from square root matrices.
#'
#' Given a list of estimates of the lower-triangular Cholesky square roots of
#' component covariance matrices, a list of true lower-triangular Cholesky
#' square roots of component covariance matrices, an estimate of the total
#' variation, and the true total variation, \code{multi_stein_loss} will
#' calculate multiway Stein's loss between the estimates and the truth.
#'
#' Multiway Stein's loss is a generalization of Stein's loss. More details on
#' multiway Stein's loss and the Bayes rules under it can be found in
#' \href{http://www.sciencedirect.com/science/article/pii/S0047259X15000330}{
#' Gerard and Hoff (2015)}.
#'
#' The function \code{multi_stien_loss_cov} also calculates multiway Stein's
#' loss, but uses the component covariance matrices (not the Cholesky roots) as
#' input.
#'
#' @param B A list of lower triangular matrices. These are the 'estimates' of
#'   the lower-triangular Cholesky square roots of the component covariance
#'   matrices.
#' @param Psi A list of lower triangular matrices. These are the 'true'
#'   lower-triangular Cholesky square roots of the component covariance
#'   matrices.
#' @param b A numeric. This is an 'estimate' of the total variation parameter,
#'   the 'standard devation' version of it.
#' @param psi A numeric. This is the 'true' total variation parameter, the
#'   'standard devation' version of it.
#'
#' @return A numeric, the multiway Stein's loss between the 'truth' and the
#'   'estimates'.
#'
#' @export
#'
#' @references Gerard, D., & Hoff, P. (2015).
#'   \href{http://www.sciencedirect.com/science/article/pii/S0047259X15000330}{Equivariant
#'   minimax dominators of the MLE in the array normal model}. \emph{Journal of
#'   Multivariate Analysis}, 137, 32-49.
#'
#' @author David Gerard.
#'
#' @keywords equivariance loss
#'
#' @seealso \code{\link{multi_stein_loss_cov}}, \code{\link{get_equi_bayes}}.
multi_stein_loss <- function(B, Psi, b, psi) {
    p <- sapply(B, nrow)
    n <- length(p)
    stein_sum <- 0
    for (index in 1:n) {
        stein_sum <- stein_sum + prod(p) / p[index] *
          tr(B[[index]] %*% t(B[[index]]) %*% solve(Psi[[index]] %*% t(Psi[[index]])))
    }
    b ^ 2 / psi ^ 2 * stein_sum - n * prod(p) * log(b ^ 2 / psi ^ 2) - n * prod(p)
}

#' Calculate multiway Stein's loss from component covariance matrices.
#'
#' Given a list of estimated component covariance matrices, a list of true
#' component covariance matrices, an estimate of the total variation, and the
#' true total variation, \code{multi_stein_loss_cov} will calculate multiway
#' Stein's loss between the estimates and the truth.
#'
#' Multiway Stein's loss is a generalization of Stein's loss. More details on
#' multiway Stein's loss and the Bayes rules under it can be found in
#' \href{http://www.sciencedirect.com/science/article/pii/S0047259X15000330}{
#' Gerard and Hoff (2015)}.
#'
#' The function \code{multi_stien_loss} also calculates multiway Stein's loss,
#' but uses the lower-triangular Cholesky square roots of the component
#' covariance matrices as input.
#'
#' @param B A list of positive definite matrices. These are the 'estimates' of
#'   the component covariance matrices.
#' @param Sigma A list of positive definite matrices. These are the 'true'
#'   component covariance matrices.
#' @param b A numeric. This is an 'estimate' of the total variation parameter,
#'   the 'standard devation' version of it.
#' @param sigma A numeric. This is the 'true' total variation parameter, the
#'   'standard devation' version of it.
#'
#' @return A numeric, the multiway Stein's loss between the 'truth' and the
#'   'estimates'.
#'
#' @export
#'
#' @references Gerard, D., & Hoff, P. (2015).
#'   \href{http://www.sciencedirect.com/science/article/pii/S0047259X15000330}{Equivariant
#'   minimax dominators of the MLE in the array normal model}. \emph{Journal of
#'   Multivariate Analysis}, 137, 32-49.
#'
#' @author David Gerard.
#'
#' @keywords equivariance loss
#'
#' @seealso \code{\link{multi_stein_loss}}, \code{\link{get_equi_bayes}}.
multi_stein_loss_cov <- function(B, Sigma, b, sigma) {
    p <- sapply(B, nrow)
    n <- length(p)
    stein_sum <- 0
    for (index in 1:n) {
        stein_sum <- stein_sum + prod(p) / p[index] * tr(B[[index]] %*% solve(Sigma[[index]]))
    }
    loss_val <- b ^ 2 / sigma ^ 2 * stein_sum - n * prod(p) * log(b ^ 2 / sigma ^ 2) - n * prod(p)
    return(loss_val)
}
