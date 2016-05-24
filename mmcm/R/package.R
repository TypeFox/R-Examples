#' The modified maximum contrast method package
#'
#' This package provides an implementation of modified maximum contrast methods and the maximum contrast method.
#' This version supports functions \code{\link{mmcm.mvt}}, \code{\link{mcm.mvt}} that gives P-value by using randomized quasi-Monte Carlo method from \code{\link[mvtnorm:pmvt]{pmvt}} function of package \code{mvtnorm}, and \code{\link{mmcm.resamp}} that gives P-value by using the permutation method.
#' In a one-way problem testing pattern of several factor level means, the maximum contrast statistics (Yoshimura, I., 1997) may be used. But under unequal sample size situations, denominator of the maximum contrast statistics is overestimated. Thus we propose a modified maximum contrast statistics for the unequal sample size situation. Denominetor of the modified maximum contrast statistics is not influenced under the unequal sample size situation.
#'
#' @name mmcm-package
#' @aliases mmcm
#' @docType package
#' @author
#' Author: Kengo NAGASHIMA and Yasunori SATO
#' 
#' Maintainer: Kengo NAGASHIMA \email{nshi@@chiba-u.jp}
#' @references
#' Nagashima, K., Sato, Y., Hamada, C. (2011).
#' A modified maximum contrast method for unequal sample sizes in pharmacogenomic studies
#' \emph{Stat Appl Genet Mol Biol.} \strong{10}(1): Article 41.
#' \url{http://dx.doi.org/10.2202/1544-6115.1560}
#' 
#' Sato, Y., Laird, N.M., Nagashima, K., et al. (2009).
#' A new statistical screening approach for finding pharmacokinetics-related genes in genome-wide studies.
#' \emph{Pharmacogenomics J.} \strong{9}(2): 137--146.
#' \url{http://www.ncbi.nlm.nih.gov/pubmed/19104505}
#' 
#' Yoshimura, I., Wakana, A., Hamada, C. (1997).
#' A performance comparison of maximum contrast methods to detect dose dependency.
#' \emph{Drug Information J.} \strong{31}: 423--432.
#' @seealso
#' \code{\link{mcm.mvt}},
#' \code{\link{mmcm.mvt}},
#' \code{\link{mmcm.resamp}}
#' @useDynLib mmcm
#' @import mvtnorm
#' @keywords htest
NULL
