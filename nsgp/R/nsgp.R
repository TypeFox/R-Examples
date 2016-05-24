#' Gaussian process regression and statistical testing
#' 
#' This package implements the non-stationary gaussian processes for 
#' one- and two-sample cases, and statistical likelihood ratio tests 
#' for distinguishing when two time-series are significantly different. 
#' The package offers two main functions: \code{\link{gpr1sample}} and
#'  \code{\link{gpr2sample}}.
#'  
#' The function \code{\link{gpr1sample}} learns a gaussian process model 
#' that uses either stationary or non-stationary gaussian kernel, which 
#' assumes a perturbation at (time) point 0. The non-stationarity is controlled
#' by a time-dependent lengthscale in the gaussian kernel. The time-dependency
#' \eqn{l - (l - l_{min})e^{-ct}}{l - (l - l_{min})e^{-ct}} follows exponential
#' decay, such that it starts at value \code{l.min} and grows logarithmically 
#' to \code{l} by curvature parameter \code{c}. 
#'
#' In \code{\link{gpr2sample}} we compare control and case time-series by 
#' building GP models for both of them individually, while also building a 
#' third null model for joint data (assume that data come from the same process). 
#' The null model and the case/control models are then compared with likelihood 
#' ratios for significant different along time. The package includes standard 
#' marginal log likelihood (MLL) ratio, and three novel ones: expected marginal log 
#' likelihood (EMLL) measures the ratio between the models while discarding data;
#' posterior concentration (PC) ratio measures the difference of variance between
#' null and individual models; and noisy posterior concentration (NPC) ratio also
#' compares observational noises.
#'
#' @docType package
#' @name nsgp
NULL

#' Toy time-series data for testing, contains two time-series, case and control
#' 
#' @source randomly generated data
#' @docType data
#' @keywords datasets
#' @name toydata
#' @usage data(toydata)
#' @format A list of two dataframes containing 24 (x,y)
#'  pairs both, i.e. toydata$ctrl and toydata$case
NULL

#' Toy time-series model for testing
#' 
#' contains the GP models for the 'toydata'
#' 
#' @name toygps
#' @docType data
#' @usage data(toygps)
#' @keywords datasets
#' @source prelearned model with \code{\link{gpr2sample}}
#' @format A \code{gppack} object
NULL

