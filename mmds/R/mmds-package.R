# No code here, just documentation!


#' A fitted Mixture Model Detection Function Object
#' 
#' The fitted mixture model detection function object returned by 
#' \code{\link{fitmix}}. Knowledge of most of this is not useful. Use 
#' \code{link{summary.ds.mixture}} for result summaries.
#'
#' @name ds.mixture
#' @section Details:
#'  A \code{\link{ds.mixture}} object has the following elements:
#'  \tabular{ll}{
#'    distance \tab Vector of distances used in the analysis.\cr
#'    likelihood \tab Value of the log-likelihood at the maxima.\cr
#'    pars \tab Parmeter estimates. See \code{\link{mmds.pars}} for more 
#'              information.\cr
#'    mix.terms \tab Number of mixture terms fit.\cr
#'    width \tab Truncation distance used.\cr
#'    z \tab List containing the matrix of covariates used. Output from 
#'            \code{\link{model.matrix}}.\cr
#'    zdim \tab Number of columns of \code{z}. See \code{\link{mmds.pars}} for 
#'              more information.\cr
#'    hessian \tab Hessian matrix at the maxima.\cr
#'    pt \tab Logical indicating whether the data were from a point transect 
#'            survey.\cr
#'    data \tab Data frame after truncation.\cr
#'    ftype \tab Type of detection function.\cr
#'    ctrl.options \tab Options passed to \code{\link{optim}}.\cr
#'    showit \tab Debug level.\cr
#'    opt.method \tab Optimisation method used.\cr
#'    usegrad \tab Were analytic gradients used?\cr
#'    model.formula \tab Model formula.\cr
#'    mu \tab Per-observation effective trip width/effective area of 
#'            detection.\cr
#'    pa.vec \tab Vector of per-observation detectabilities.\cr
#'    N \tab Estimate of N in the covered area (Horvitz-Thompson).\cr
#'    pa \tab Average detectability.\cr
#'    pars.se \tab Standard errors of the parameters.\cr
#'    N.se \tab Standard error of the Horvitz-Thompson estimate of the 
#'              abundance.\cr
#'    pa.se \tab Standard error of the average detectability.\cr
#'    aic \tab AIC of the fitted model.\cr
#'    cvm \tab Cramer-von Mises GoF test results. List containing: \code{p}, 
#'              the p-value and \code{W}, the test statistic.\cr
#'    ks \tab Kolmogorov-Smirnov test results. List containing: \code{p}, the 
#'            p-value and \code{Dn}, the test statistic. See 
#'            \code{\link{mmds.gof}} for more information.\cr}
#'
#' @author David L. Miller
#' @section Note:
#' \code{ds.mixture} objects can be passed to \code{\link{step.ds.mixture}} to 
#'    select number of mixture components based on AIC score. 
NULL

#' Goodness of fit for mixture model detection functions
#'
#' Goodness of fit testing for detection for mixture model detection functions.
#' 
#' @name mmds.gof
#' @section Details:
#' Two goodness of fit tests are provided: the Cramer-von Mises and the 
#' Kolmogorov-Smirnov. Both are implemented as in Buckland et al. (2004).
#' 
#' Print methods are provided, so accessing the \code{ks} and \code{cvm} 
#' elements of a \code{\link{ds.mixture}} object will give suitable summaries.
#'
#' David L. Miller
NULL

#' Parameters in \code{mmds}
#' 
#' The internal parametrisation used in \code{mmds} is not directly 
#' interpretable. This man page aims to explain how to interpret the parameters
#' and transform them into useful information.
#'
#' @name mmds.pars
#' @section Details:
#' Parametrisation works differently for the scale parameters of the 
#' half-normals and for the mixture proportions.
#' 
#' The scale parameters of the half-normals (or their constituent parameters in
#' the case of a covariate model) are given on the log scale.
#' 
#' The mixture proportions are transformed to a parametrisation that allows 
#' values to lie over all of the real line (see Miller and Thomas for details).
#' 
#' The parameter vector is made up of the scale parameters followed by the 
#' mixture parameters. In the non-covariate case the former is the length of 
#' the number of mixtures (\code{mix.terms}) and the latter is of length 
#' \code{mix.terms}-1. When the model has covariates the scale paramters are 
#' given as \code{mix.terms} intercepts followed by the covariate parameters.
#' 
#' The function \code{\link{getpars}} transforms the parameters (\code{$pars} 
#' element) in a \code{\link{ds.mixture}} object to a named list.
#' 
#' Calling \code{\link{summary.ds.mixture}} will show the mixture proportions.
#' 
#' The parameter \code{initialvalues} supplied to \code{\link{fitmix}}.
#' 
#' @author David L. Miller
NULL
