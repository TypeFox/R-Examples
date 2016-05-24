#' Mandallaz' model-assisted small area estimators
#' 
#' an S4 implementation of the unbiased extension of the 
#' model-assisted' synthetic-regression estimator proposed by 
#' Mandallaz (2013), Mandallaz et al. (2013) and Mandallaz (2014).\cr
#' It yields smaller variances than the standard bias correction, 
#' the generalised regression estimator.
#' 
#' This package provides Mandallaz' extended synthetic-regression estimator for two- and
#' three-phase sampling designs with or without clustering.\cr
#' See vignette('maSAE', package = 'maSAE') and demo('maSAE', package = 'maSAE') for
#' introductions, \code{"\link[=saeObj-class]{class?maSAE::saeObj}"} and
#' \code{"\link[=predict]{?maSAE::predict}"} for help on the main feature.
#' 
#' @note Model-assisted estimators use models to improve the efficiency (i.e. reduce
#' prediction error compared to design-based estimators) but need not assume them to be
#' correct as in the model-based approach, which is advantageous in official
#' statistics.  
#' @name maSAE-package
#' @aliases maSAE-package tmp.maSAE.skeleton
#' @docType package
#' @seealso There are a couple packages for model-\strong{based} small area estimation, see 
#' \code{\link[sae:sae-package]{sae}}, 
#' \code{\link[rsae:rsae-package]{rsae}}, 
#' hbsae and
#' \code{\link[JoSAE:JoSAE-package]{JoSAE}}.
#' @references 
#' \cite{
#' Mandallaz, D. 2013 
#' Design-based properties of some small-area estimators in forest
#' inventory with two-phase sampling. 
#' Canadian Journal of Forest Research \bold{43}(5), pp. 441--449. 
#' doi: \href{http://dx.doi.org/10.1139/cjfr-2012-0381}{10.1139/cjfr-2012-0381}.
#' }
#'
#' \cite{
#' Mandallaz, and Breschan, J.  and  Hill, A. 2013 
#' New regression estimators in forest inventories with two-phase sampling and partially 
#' exhaustive information: a design-based Monte Carlo approach with applications to 
#' small-area estimation. 
#' Canadian Journal of Forest Research \bold{43}(11), pp. 1023--1031. 
#' doi: \href{http://dx.doi.org/10.1139/cjfr-2013-0181}{10.1139/cjfr-2013-0181}.
#' }
#'
#' \cite{
#' Mandallaz, D. 2014 
#' A three-phase sampling extension of the generalized regression
#' estimator with partially exhaustive information. 
#' Canadian Journal of Forest Research \bold{44}(4), pp. 383--388. 
#' doi: \href{http://dx.doi.org/10.1139/cjfr-2013-0449}{10.1139/cjfr-2013-0449}.
#' }
#' 
#' @keywords package
#' @examples
#' \dontrun{vignette('maSAE', package = 'maSAE')}
#' \dontrun{demo('design', package = 'maSAE')}
#' \dontrun{demo('maSAE', package = 'maSAE')}
#' 
NULL
