#' \code{Asia} dataset.
#' 
#' The \code{Asia} dataset contains 10000 complete (no missing data, no latent variables) randomly generated items of the \code{Asia} Bayesian Network.
#' No imputation needs to be performed, so only raw data is present.
#' 
#' The data the BNDataset object is built from is located in files \code{pkg_folder/extdata/asia_10000.header} and \code{pkg_folder/extdata/asia_10000.data}.
#'
#' @name asia_10000
#' @rdname asia_10000
#' @docType data
#' 
#' @references
#' S. Lauritzen, D. Spiegelhalter. Local Computation with Probabilities on Graphical Structures and their Application to Expert Systems
#' (with discussion). Journal of the Royal Statistical Society: Series B (Statistical Methodology), 50(2):157-224, 1988.
#' 
#' @format a \code{\link{BNDataset}} with raw data slow filled.
#' 
#' @seealso \code{\link{asia}}
NULL
# "asia_10000"

#' \code{Child} dataset.
#' 
#' The \code{Child} dataset contains 5000 randomly generated items with missing data (no latent variables) of the \code{Child} Bayesian Network.
#' Imputation is performed, so both raw and imputed data is present.
#'
#' The data the BNDataset object is built from is located in files \code{pkg_folder/extdata/extdata/Child_data_na_5000.header} and \code{pkg_folder/extdata/extdata/Child_data_na_5000.data}.
#'
#' @name child_NA_5000
#' @rdname child_NA_5000
#' @docType data
#' 
#' @references
#' D. J. Spiegelhalter, R. G. Cowell (1992). Learning in probabilistic expert systems. In Bayesian Statistics 4
#' (J. M. Bernardo, J. 0. Berger, A. P. Dawid and A. F. M. Smith, eds.) 447-466. Clarendon Press, Oxford. 
#' 
#' @format a \code{\link{BNDataset}} with a raw and imputed data slow filled with 5000 items.
#' 
#' @seealso \code{\link{child}}
NULL
# "child_NA_5000"
