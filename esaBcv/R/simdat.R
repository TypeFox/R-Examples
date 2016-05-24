#' Example Dataset
#'
#' The data is a simulated data set where the data matrix is generated from the latent 
#' factor model
#' \deqn{Y = n^{1/2}U D V' + E \Sigma^{1/2}}
#' where \eqn{D} and \eqn{\Sigma} are diagonal matrices, and \eqn{U} and \eqn{V} 
#' are orthogonal. \eqn{V'} means _V transposed_. For the factors, we include one giant
#' factor, five useful factors, one harmful factor and one undetectable factor.
#' For more details of the simulation method
#' used, please refer to Appendix A.1 of Owen and Wang (2015) Bi-cross-validation for factor analysis, \url{http://arxiv.org/abs/1503.03515}.
#'
#' The dataset is a list of components:
#' \itemize{
#' \item\code{Y} a data matrix of 200 by 1000, where each row is a sample and each column 
#' is a variable
#' \item\code{U} the orthogonal factor matrix \eqn{U} of size 200 by 8.
#' \item\code{V} the orthogonal factor matrix \eqn{V} of size 1000 by 8.
#' \item\code{D} the vector of diagonal entries of \eqn{D}.
#' \item\code{Sigma} the vector of diagonal entries of \eqn{\Sigma}.
#' \item\code{oracle.r} the oracle rank (the optimal number of factors that should be kept)
#' of the factor matrix.
#' }
#'
#' @name simdat
#' @docType data
NULL
