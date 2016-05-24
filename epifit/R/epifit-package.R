#' Flexible Modelling Functions for Epidemiological Data Analysis.
#'
#' Provides flexible model fitting used in epidemiological data analysis
#' by a unified model specification, along with some data manipulation functions.
#' This package covers fitting of variety models including Cox regression models,
#' linear regression models, Poisson regression models, logistic models and
#' others whose likelihood is expressed in negative binomial, gamma and Weibull
#' distributions.
#' 
#' @name epifit-package
#' @aliases epifit-pkg
#' @docType package
#' @author
#' Author: Kazutaka Doi, Kei Sakabe and Masataka Taguri
#' Maintainer: Kazutaka Doi \email{kztkdi@@gmail.com}
#' @seealso
#' \code{\link{AIC.epifit}},
#' \code{\link{calcAge}},
#' \code{\link{convertFromFactor}},
#' \code{\link{convertNA}},
#' \code{\link{countNA}},
#' \code{\link{epifit}},
#' \code{\link{extractVariable}},
#' \code{\link{listNumericIncompatibility}},
#' \code{\link{modes}},
#' \code{\link{pullOneValue}},
#' \code{\link{pytable}},
#' \code{\link{removeVariable}},
#' \code{\link{showContents}}
#' @importFrom stats dbinom dgamma dnbinom dnorm dpois dweibull integrate na.exclude na.fail na.omit na.pass nlm optim pchisq pnorm
#' @importFrom MASS ginv
#' @keywords models
NULL
