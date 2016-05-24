#' Functions for Extreme Value Mixture Modelling, Threshold Estimation and
#' Boundary Corrected Kernel Density Estimation
#'
#' \tabular{ll}{
#' Package: \tab evmix\cr
#' Type: \tab Package\cr
#' Version: \tab 0.2-6\cr
#' Date: \tab 2015-05-27\cr
#' License: \tab GPL-3\cr
#' LazyLoad: \tab yes\cr
#' }
#'
#' The usual distribution functions, maximum likelihood inference and model
#' diagnostics for univariate stationary extreme value mixture models are
#' provided.
#' 
#' Kernel density estimation including various boundary corrected kernel density
#' estimation methods and a wide choice of kernels, with cross-validation
#' likelihood based bandwidth estimators are included.
#' 
#' Reasonable consistency with the base functions in the \code{evd} package is
#' provided, so that users can safely interchange most code.
#'  
#' @name evmix-package
#' @aliases evmix
#' @docType package
#' @title Extreme Value Mixture Modelling, Threshold Estimation and Boundary Corrected Kernel Density Estimation
#' @author Carl Scarrott and Yang Hu, University of Canterbury, New Zealand \email{carl.scarrott@@canterbury.ac.nz}
#' @references
#' \url{http://www.math.canterbury.ac.nz/~c.scarrott/evmix}
#' 
#' Scarrott, C.J. and MacDonald, A. (2012). A review of extreme value
#' threshold estimation and uncertainty quantification. REVSTAT - Statistical
#' Journal 10(1), 33-59. Available from \url{http://www.ine.pt/revstat/pdf/rs120102.pdf}
#' 
#' Hu, Y. (2013). Extreme value mixture modelling: An R package and simulation study.
#' MSc (Hons) thesis, University of Canterbury, New Zealand.
#' \url{http://ir.canterbury.ac.nz/simple-search?query=extreme&submit=Go}
#' 
#' MacDonald, A. (2012). Extreme value mixture modelling with medical and
#' industrial applications. PhD thesis, University of Canterbury, New Zealand.
#' \url{http://ir.canterbury.ac.nz/bitstream/10092/6679/1/thesis_fulltext.pdf}
#' 
#' @seealso \code{\link[evd:gpd]{evd}}, \code{\link[ismev:ismev]{ismev}} and 
#' \code{\link[condmixt:condmixt-package]{condmixt}}
#' @import stats graphics MASS splines gsl
#' @importFrom SparseM as.matrix.csr rbind.matrix.csr slm.wfit is.matrix.csr
NULL
