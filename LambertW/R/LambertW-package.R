#' @title R package for Lambert W\eqn{ \times} F distributions
#' @name LambertW-package
#' @aliases LambertW
#' @docType package
#' 
#' @description
#'
#' This package is based on notation, definitions, and results of Goerg (2011,
#'     2015, 2016).  I will not include these references in the description of
#'     each single function.
#' 
#' Lambert W\eqn{ \times} F distributions are a general framework to model and
#'     transform skewed, heavy-tailed data. Lambert W\eqn{ \times} F random
#'     variables (RV) are based on an input/ouput system with input RV \eqn{X
#'     \sim F_X(x \mid \boldsymbol \beta)} and output \eqn{Y}, which is a
#'     non-linearly transformed version of X -- with similar properties to X,
#'     but slightly skewed and/or heavy-tailed.  Then Y has a 'Lambert W
#'     \eqn{\times F_X}' distribution - see References.
#' 
#' \code{\link{get_distnames}} lists all implemented Lambert W \eqn{\times} F
#'     distributions in this package. If you want to generate a
#'     skewed/heavy-tailed version of a distribution that is not implemented,
#'     you can use the do-it-yourself modular toolkit
#'     (\code{\link{create_LambertW_input}} and
#'     \code{\link{create_LambertW_output}}). It allows users to quickly
#'     implement their own Lambert W x 'MyFavoriteDistribution' and use it in
#'     their analysis right away.
#' 
#' This package contains several functions to analyze skewed and heavy-tailed
#'     data: simulate random samples (\code{\link{rLambertW}}), evaluate pdf and
#'     cdf (\code{\link{dLambertW}} and \code{\link{pLambertW}}), estimate
#'     parameters (\code{\link{IGMM}} and \code{\link{MLE_LambertW}}), compute
#'     quantiles (\code{\link{qLambertW}}), and plot/print results nicely
#'     (\code{\link{plot.LambertW_fit}}, \code{\link{print.LambertW_fit}},
#'     \code{\link{summary.LambertW_fit}}).
#' 
#' Probably the most useful function is \code{\link{Gaussianize}}, which works
#'     similarly to \code{\link[base]{scale}}, but makes your data Gaussian (not
#'     just centers and scales it, but also makes it symmetric and removes
#'     excess kurtosis).
#' 
#' If you use this package in your work please cite it
#'     (\code{citation("LambertW")}).  You can also send me an implementation of
#'     your 'Lambert W \eqn{\times} YourFavoriteDistribution' to add to the
#'     \pkg{LambertW} package (and I will reference your work introducing your
#'     'Lambert W \eqn{\times} YourFavoriteDistribution' here.)
#' 
#' Feel free to contact me for comments, suggestions, code improvements,
#'     implementation of new input distributions, bug reports, etc.
#' 
#' @author Author and maintainer: Georg M. Goerg (im (at) gmge.org)
#' @references
#' Goerg, G.M. (2011). \dQuote{Lambert W Random Variables - A New Family of
#'     Generalized Skewed Distributions with Applications to Risk
#'     Estimation}. Annals of Applied Statistics, 5 (3), 2197-2230.
#'     (\url{http://arxiv.org/abs/0912.4554}).
#' 
#' Goerg, G.M. (2015). \dQuote{The Lambert Way to Gaussianize heavy-tailed data
#'     with the inverse of Tukey's h transformation as a special case}.  The
#'     Scientific World Journal: Probability and Statistics with Applications in
#'     Finance and Economics. Available at
#'     \url{http://www.hindawi.com/journals/tswj/aa/909231/}.
#'
#' Goerg, G.M. (2016).  \dQuote{Rebuttal of the ``Letter to the Editor of
#'     Annals of Applied Statistics'' on Lambert W x F distributions and the
#'     IGMM algorithm}.  Available on arxiv.
#' 
#' @keywords package
#' @import MASS stats graphics
#' @useDynLib LambertW
#' @importFrom Rcpp sourceCpp evalCpp
#' @examples
#'  
#' # Replicate parts of the analysis in Goerg (2011)
#' data(AA)
#' y <- AA[AA$sex=="f", "bmi"]
#' test_normality(y)
#' 
#' fit.gmm <- IGMM(y, type = "s")
#' summary(fit.gmm)  # gamma is significant and positive
#' plot(fit.gmm)
#' 
#' # Compare empirical to theoretical moments (given parameter estimates)
#' moments.theory <- 
#'  mLambertW(theta = list(beta = fit.gmm$tau[c("mu_x", "sigma_x")], 
#'                         gamma = fit.gmm$tau["gamma"]), 
#'            distname = "normal")
#' TAB <- rbind(unlist(moments.theory),
#'              c(mean(y), sd(y), skewness(y), kurtosis(y)))
#' rownames(TAB) <- c("Theoretical (IGMM)", "Empirical")
#' TAB
#' 
#' x <- get_input(y, fit.gmm$tau)
#' test_normality(x) # input is normal -> fit a Lambert W x Gaussian by MLE
#' 
#' fit.ml <- MLE_LambertW(y, type = "s", distname = "normal", hessian = TRUE)
#' summary(fit.ml)
#' plot(fit.ml)
#'
NULL
