#' Random Effects and Sample Selection Models for Panel Counting Data
#' @description  A high performance package for estimating counting models with random effects and sample selection in panel counting data, namely counting data with repeated observations on individuals over time. \cr
#' @section Functions:
#' ProbitRE: Probit model with random effects on individuals \cr \cr
#' PoissonRE: Poisson model with random effects on individuals \cr \cr
#' PLN_RE: Poisson Lognormal model with random effects on individuals \cr \cr
#' CRE: PoissonRE and ProbitRE model with correlated random effects on individuals \cr \cr
#' CRE_SS: PLN_RE and ProbitRE model with correlated random effects on individual level and correlated error terms on <individual, time> level \cr \cr
#' @docType package
#' @name PanelCount
#' @importFrom statmod gauss.quad
#' @importFrom Rcpp evalCpp
#' @importFrom stats binomial dnorm dpois glm model.frame model.matrix model.response optim pchisq pnorm poisson runif
#' @importFrom utils tail
#' @useDynLib PanelCount
#' @references 1. Jing Peng and Christophe Van den Bulte. Participation vs. Effectiveness of Paid Endorsers in Social Advertising Campaigns: A Field Experiment. Working Paper.
#' @references 2. Jing Peng and Christophe Van den Bulte. How to Better Target and Incent Paid Endorsers in Social Advertising Campaigns: A Field Experiment. In Proceedings of the 2015 International Conference on Information Systems.
NULL
#> NULL