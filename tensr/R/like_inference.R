#' Calculate the likelihood ratio test statistic.
#'
#' Calulate the likelihood ratio test statistic for Kronecker structured
#' covariance models.
#'
#' The LRT statistic is the exact same for all elliptically distributed
#' Kronecker structured covariance models (not just the normal). The
#' distribution of the likelihood ratio test statistic does change.
#'
#' @param sig_null A numeric. The MLE of the total variation parameter under the
#'   null (the standard deviation version).
#' @param sig_alt A numeric. The MLE of the total variation parameter under the
#'   alternative (the standard deviation version).
#' @param p A vector of integers. The dimension of the array.
#'
#' @return A numeric. The likelihood ratio test statistic.
#'
#' @author David Gerard.
#'
#' @keywords likelihood
#'
#' @export
#'
#' @references Gerard, D. C., & Hoff, P. D. (2014).
#'   \href{http://arxiv.org/abs/1410.1094}{A higher-order LQ decomposition for
#'   separable covariance models}. \emph{arXiv preprint arXiv:1410.1094.}
#'
#' @seealso \code{\link{holq}} for obtaining the MLE of the total variation
#'   parameter.
#'
#'   \code{\link{lrt_null_dist_dim_same}} for getting the null distribution of
#'   the likelihood ratio test statistic.
#'
lrt_stat <- function(sig_null, sig_alt, p) {
    return((log(sig_null) - log(sig_alt)) * prod(p) * 2)
}

#' Calculate the AIC and BIC.
#'
#' Calculate the AIC and BIC for Kronecker structured covariance models,
#' assuming the array normal distribution.
#'
#' The AIC and BIC depend only on the data through the MLE of the total
#' variation parameter. Given this, the dimension of the array, and a
#' specification of which modes are the identity and which are unstructured,
#' this function will calculate the AIC and BIC.
#'
#' @param sig_squared A numeric. The MLE of sigma^2 in the array normal model
#'   (the 'variance' form of the total variation parameter).
#' @param p A vector of integers. The dimension of the data array (including
#'   replication modes).
#' @param mode_ident A vector of integers. The modes assumed to have identity
#'   covariances.
#' @param mode_diag A vector of integers. The modes assumed to have diagional
#'   covariances.
#' @param mode_unstructured A vector of integers. The modes of assumed to have
#'   unstructured covariances.
#'
#' @export
#'
#' @seealso \code{\link{holq}} for obtaining \code{sig_squared}.
#'
#' @author David Gerard.
#'
#' @keywords likelihood
#'
#' @return \code{AIC} A numeric. The AIC of the model.
#'
#'   \code{BIC} A numeric. The BIC of the model.
#'
#' @examples
#' # Generate random array data with first mode having unstructured covariance
#' #  second having diagonal covariance structure and third mode having identity
#' #  covariance structure.
#' set.seed(857)
#' p <- c(4, 4, 4)
#' Z <- array(stats::rnorm(prod(p)), dim = p)
#' Y <- atrans(Z, list(tensr:::rwish(diag(p[1])), diag(1:p[2]), diag(p[3])))
#'
#' # Use holq() to fit various models.
#' false_fit1 <- holq(Y, mode_rep = 1:3) ## identity for all modes
#' false_fit2 <- holq(Y, mode_rep = 2:3) ## unstructured first mode
#' true_fit <- holq(Y, mode_rep = 3, mode_diag = 2) ## correct model
#'
#' # Get AIC and BIC values.
#' false_aic1 <- array_bic_aic(false_fit1$sig ^ 2, p, mode_ident = 1:length(p))
#' false_aic2 <- array_bic_aic(false_fit2$sig ^ 2, p, mode_ident = 2:length(p),
#'                             mode_unstructured = 1)
#' true_aic <- array_bic_aic(true_fit$sig ^ 2, p, mode_ident = 2:length(p), mode_diag = 1)
#'
#' # Plot the results.
#' plot(c(false_aic1$AIC, false_aic2$AIC, true_aic$AIC), type = "l",
#'      xaxt = "n", xlab = "Model", ylab = "AIC", main = "AIC")
#' axis(side = 1, at = 1:3, labels = c("Wrong Model 1", "Wrong Model 2", "Right Model"))
#'
#' plot(c(false_aic1$BIC, false_aic2$BIC, true_aic$BIC), type = "l", xaxt = "n",
#'      xlab = "Model", ylab = "BIC", main = "BIC")
#' axis(side = 1, at = 1:3, labels = c("Wrong Model 1", "Wrong Model 2", "Right Model"))
array_bic_aic <- function(sig_squared, p, mode_ident = NULL, mode_diag = NULL, mode_unstructured = NULL) {

    if (!is.null(mode_ident)) {
        n <- prod(p[mode_ident])
    } else {
        n <- 1
    }

    if (!is.null(mode_diag)) {
        k1 <- sum(p[mode_diag] - 1)
    } else {
        k1 <- 0
    }

    if (!is.null(mode_unstructured)) {
        k2 <- sum(p[mode_unstructured] * (p[mode_unstructured] + 1) / 2 - 1)
    } else {
        k2 <- 0
    }

    k <- k1 + k2

    minus_2_log_L_hat <- prod(p) * (log(2 * pi) + log(sig_squared) + 1)
    bic <- minus_2_log_L_hat + k * log(n)
    aic <- minus_2_log_L_hat + 2 * k
    return(list(AIC = aic, BIC = bic))
}

#' Draw from null distribution of likelihood ratio test statistic.
#'
#' When testing for the covariance structure of modes, this function may be used
#' to draw a sample from the null distribution of the likelihood ratio test
#' stistics, whose distribution doesn't depend on any unknown parameters under
#' the null.
#'
#' Let \eqn{vec(X)} be \eqn{N(0,\Sigma)}. Given two nested hypotheses, \deqn{H_1:
#' \Sigma = \Psi_K\otimes\cdots\otimes\Psi_1} versus \deqn{H_0: \Sigma =
#' \Omega_K\otimes\cdots\otimes\Omega_1,} this function will draw from the null
#' distribution of the likelihood ratio test statistic. The possible options are
#' that \eqn{\Psi_i} or \eqn{\Omega_i} are the identity matrix, a diagonal
#' matrix, or any positive definite matrix. By default, it's assumed that these
#' matrices are any positive definite matrix.
#'
#' Unfortunately, this fuction does not support testing for the hypothesis of
#' modeling the covariance between two modes with a single covariance matrix. I
#' might code this up in later versions.
#'
#' @param p A vector of integers. The dimensions of the array.
#' @param null_ident A vector of integers. The modes that under the null have
#'   identity covariance.
#' @param alt_ident A vector of integers. The modes that under the alternative
#'   have the identity covariance.
#' @param null_diag A vector of integers. The modes that under the null have
#'   diagonal covariance.
#' @param alt_diag A vector of integers. The modes that under the alternative
#'   have diagonal covariance.
#' @param reference_dist Two options are supported, 'normal' and 't'. If 't' is
#'   specified, you have to specify \code{t_df}.
#' @param t_df A numeric. If \code{reference_dist} is 't', then this is the
#'   degrees of freedom of the t_distribution that the array is distributed
#'   under.
#' @param itermax An integer. The number of draws from the null distribution of
#'   the likelihood ratio test statistic that is to be performed.
#' @param holq_itermax An integer. The maximum number of block coordinate ascent
#'   iterations to perform when calculating the MLE at each step.
#' @param holq_tol A numeric. The stopping criterion when calculating the MLE.
#'
#' @return A vector of draws from the null distribution of the likelihood ratio
#'   test statistic.
#'
#' @seealso \code{\link{lrt_stat}} for calculating the likelihood ratio test
#'   statistic.
#'
#' @export
#'
#' @keywords likelihood
#'
#' @references Gerard, D. C., & Hoff, P. D. (2014).
#'   \href{http://arxiv.org/abs/1410.1094}{A higher-order LQ decomposition for
#'   separable covariance models}. \emph{arXiv preprint arXiv:1410.1094.}
#'
#' @author David Gerard.
#'
#' @examples
#' #Test for all identity versus all unconstrained.
#' p = c(4,4,4)
#' null1 <- lrt_null_dist_dim_same(p,null_ident = 1:3)
#'
#' #Generate Null Data
#' X <- array(stats::rnorm(prod(p)), dim = p)
#' sig_null <- holq(X, mode_rep = 1:3)$sig
#' sig_alt <- holq(X)$sig
#' lrt_x <- lrt_stat(sig_null, sig_alt, p = p)
#' p_value <- mean(null1 > lrt_x)
#'
#' hist(null1,main = 'Null Distribution of LRT', xlab = 'LRT Statistic')
#' abline(v = lrt_x, lty = 2, col = 2, lwd = 2)
#' legend('topleft', 'Observed LRT Statistic', lty = 2, col = 2, lwd = 2)
#' mtext(side = 1, paste('P-value = ', round(p_value, digits = 2), sep = ''),
#'       line = 2)
#'
#' #-------------------------------------------------------------------------
#'
#' #Test for all identity versus all mode 1 identity,
#' #  mode 2 diagonal, mode 3 unconstrained.
#' p = c(4,4,4)
#' null2 <- lrt_null_dist_dim_same(p,null_ident = 1:3,
#'                                 alt_ident = 1, alt_diag = 2)
#'
#' #Generate Null Data
#' X <- array(stats::rnorm(prod(p)), dim = p)
#' sig_null <- holq(X, mode_rep = 1:3)$sig
#' sig_alt <- holq(X, mode_rep = 1, mode_diag = 2)$sig
#' lrt_x <- lrt_stat(sig_null, sig_alt, p = p)
#' p_value <- mean(null2 > lrt_x)
#'
#' hist(null2,main = 'Null Distribution of LRT', xlab = 'LRT Statistic')
#' abline(v = lrt_x, lty = 2, col = 2, lwd = 2)
#' legend('topleft', 'Observed LRT Statistic', lty = 2, col = 2, lwd = 2)
#' mtext(side = 1, paste('P-value = ', round(p_value, digits = 2), sep = ''),
#'       line = 2)
lrt_null_dist_dim_same <- function(p, null_ident = NULL, alt_ident = NULL, null_diag = NULL,
                                   alt_diag = NULL, reference_dist = "normal", t_df = NULL,
                                   itermax = 100, holq_itermax = 100, holq_tol = 10 ^ -9) {
    lrt_null_vec <- rep(NA, length = itermax)
    for (index in 1:itermax) {
        if (reference_dist == "normal") {
            X <- array(stats::rnorm(prod(p)), dim = p)
        } else if (reference_dist == "t") {
            X <- array(stats::rt(prod(p), df = t_df), dim = p)
        }

        sig_null <- holq(X, mode_rep = null_ident, mode_diag = null_diag, print_diff = FALSE,
                         itermax = holq_itermax, tol = holq_tol)$sig
        sig_alt <- holq(X, mode_rep = alt_ident, mode_diag = alt_diag, print_diff = FALSE,
                        itermax = holq_itermax, tol = holq_tol)$sig
        lrt_null_vec[index] <- lrt_stat(sig_null = sig_null, sig_alt = sig_alt, p = p)
        cat(lrt_null_vec[index], "\n")
    }
    return(lrt_null_vec)
}

#' Get MLE from output of \code{holq}.
#'
#' From the output of \code{holq}, this function will calculate the
#' MLEs for the component covariance matrices and for the total
#' variation parameter.
#'
#' The function simply takes the \code{A[[i]]} output of \code{holq}
#' and returs \code{A[[i]] \%*\% t(A[[i]])}. The estimate of the total
#' variation parameter is \code{sqrt(sig ^ 2 / prod{p})}, whre \code{p} is the
#' vector of dimensions of the data array and \code{sig} is the output
#' from \code{holq}.
#'
#' @param holq_obj The output returned from \code{holq}.
#'
#' @return \code{cov_mle} A list of positive definite matrices. These
#'     are the MLEs for the component covariance matrices.
#'
#'   \code{sig_mle} A numeric. This is an estimate of the "standard
#'   deviation" form of the total variation parameter.
#'
#' @export
#'
#' @keywords likelihood
#'
#' @author David Gerard.
#'
#' @references Gerard, D. C., & Hoff, P. D. (2014).
#'   \href{http://arxiv.org/abs/1410.1094}{A higher-order LQ decomposition for
#'   separable covariance models}. \emph{arXiv preprint arXiv:1410.1094.}
#'
#' @seealso \code{\link{holq}}.
mle_from_holq <- function(holq_obj){
  p <- dim(holq_obj$Z)
  cov_mle <- listprod(holq_obj$A, lapply(holq_obj$A, t))
  sig_mle <- sqrt(holq_obj$sig ^ 2 / prod(p))  #Estimate of standard deviation, not variance.
  return(list(cov_mle = cov_mle, sig_mle = sig_mle))
}
