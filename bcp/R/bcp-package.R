#' Bayesian Analysis of Change Point Problems
#' 
#' Provides an implementation of the Barry and Hartigan (1993) product partition model for the normal errors change point problem using Markov Chain Monte Carlo.  It also (i) extends the methodology to regression models on a connected graph (Wang and Emerson, 2015) and (ii) allows estimation of change point models with multivariate responses. Parallel MCMC, previously available in bcp v.3.0.0, is currently not implemented.
#' @docType package
#' @name bcp-package
#' @author Xiaofei Wang, Chandra Erdman, and John W. Emerson
#' @import methods grid  
#' @importFrom Rcpp evalCpp
#' @useDynLib bcp
NULL
