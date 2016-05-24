#' tensr: A package for Kronecker structured covariance inference.
#'
#' This package provides a collection of functions for likelihood and
#' equivariant inference for covariance matrices under the array
#' normal model.  Also included are functions for calculating tensor
#' decompositions that are related to likelihood inference in the
#' array normal model.
#'
#' @section Introduction: Let \eqn{X} be a multidimensional array
#'     (also called a tensor) of K dimensions. This package provides a
#'     series of functions to perform statistical inference when
#'     \deqn{vec(X) \sim N(0,\Sigma),} where \eqn{\Sigma} is assumed to
#'     be Kronecker structured. That is, \eqn{\Sigma} is the Kronecker
#'     product of \eqn{K} covariance matrices, each of which has the
#'     interpretation of being the covariance of \eqn{X} along its
#'     \eqn{k}th mode, or dimension.
#'
#'   Pay particular attention to the zero mean assumption. That is,
#'   you need to de-mean your data prior to applying these
#'   functions. If you have more than one sample, \eqn{X_i} for \eqn{i
#'   = 1,\ldots,n}, then you can concatenate these tensors along a
#'   \eqn{(K+1)}th mode to form a new tensor \eqn{Y} and apply the
#'   \code{demean_tensor()} function to Y which will return a tensor
#'   that satisfies the mean-zero assumption.
#'
#'   The details of the methods in this package can be found in
#'   \href{http://www.sciencedirect.com/science/article/pii/S0047259X15000330}{Gerard
#'   and Hoff (2015)} and \href{http://arxiv.org/abs/1410.1094}{Gerard
#'   and Hoff (2014)}.
#'
#' @section Tensr functions: \code{\link{amprod}} \eqn{k}-mode product.
#'
#'   \code{\link{anorm_cd}} Array normal conditional distributions.
#'
#'   \code{\link{array_bic_aic}} Calculate the AIC and BIC.
#'
#'   \code{\link{arrIndices}} Array indices.
#'
#'   \code{\link{atrans}} Tucker product.
#'
#'   \code{\link{collapse_mode}} Collapse multiple modes into one mode.
#'
#'   \code{\link{convert_cov}} Convert the output from \code{equi_mcmc} to
#'   component covariance matrices.
#'
#'   \code{\link{demean_tensor}} Demeans array data.
#'
#'   \code{\link{equi_mcmc}} Gibbs sampler using an invariant prior.
#'
#'   \code{\link{fnorm}} Frobenius norm of an array.
#'
#'   \code{\link{get_equi_bayes}} Get the Bayes rule under multiway Stein's
#'   loss.
#'
#'   \code{\link{get_isvd}} Calculate the incredible SVD (ISVD).
#'
#'   \code{\link{holq}} Calculate the incredible higher-order LQ decomposition
#'   (HOLQ).
#'
#'   \code{\link{hooi}} Calculate the higher-order orthogonal iteration (HOOI).
#'
#'   \code{\link{hosvd}} Calculate the (truncated) higher-order SVD (HOSVD).
#'
#'   \code{\link{Kom}} Commutation matrix.
#'
#'   \code{\link{ihop}} The incredible higher-order polar decomposition (IHOP).
#'
#'   \code{\link{ldan}} Log-likelihood of array normal model.
#'
#'   \code{\link{listprod}} Element-wise matrix products between two lists.
#'
#'   \code{\link{lq}} LQ decomposition.
#'
#'   \code{\link{lrt_null_dist_dim_same}} Draw from null distribution of
#'   likelihood ratio test statistic.
#'
#'   \code{\link{lrt_stat}} Calculate the likelihood ratio test statistic.
#'
#'   \code{\link{mat}} Unfold a matrix.
#'
#'   \code{\link{mhalf}} The symmetric square root of a positive definite
#'   matrix.
#'
#'   \code{\link{mle_from_holq}} Get MLE from output of \code{holq}.
#'
#'   \code{\link{multi_stein_loss}} Calculate multiway Stein's loss from square
#'   root matrices.
#'
#'   \code{\link{multi_stein_loss_cov}} Calculate multiway Stein's loss from
#'   component covariance matrices.
#'
#'   \code{\link{multiway_takemura}} Calculate a truncated multiway Takemura
#'   estimator.
#'
#'   \code{\link{polar}} The left polar decomposition.
#'
#'   \code{\link{qr2}} QR Decomposition.
#'
#'   \code{\link{random_ortho}} Generate a list of orthogonal matrices drawn
#'   from Haar distribution.
#'
#'   \code{\link{rmirror_wishart}} Sample from the mirror-Wishart distribution.
#'
#'   \code{\link{sample_sig}} Update for total variation parameter in
#'   \code{equi_mcmc}.
#'
#'   \code{\link{sample_right_wishart}}  Gibbs update of \code{Phi_inv}.
#'
#'   \code{\link{start_ident}} Get list of identity matrices.
#'
#'   \code{\link{start_resids}} Sample covariance matrices for each mode.
#'
#'   \code{\link{tsum}} Tucker sum.
#'
#'   \code{\link{tr}} Trace of a matrix.
#'
#'   \code{\link{trim}} Truncates small numbers to 0.
#'
#' @references Gerard, D., & Hoff, P. (2015).
#'     \href{http://www.sciencedirect.com/science/article/pii/S0047259X15000330}{Equivariant
#'     minimax dominators of the MLE in the array normal
#'     model}. \emph{Journal of Multivariate Analysis}, 137, 32-49.
#'
#'     Gerard, D. C., & Hoff, P. D. (2014).
#'     \href{http://arxiv.org/abs/1410.1094}{A higher-order LQ
#'     decomposition for separable covariance models}. \emph{arXiv
#'     preprint arXiv:1410.1094.}
#'
#' @docType package
#' @name tensr
NULL
# > NULL
