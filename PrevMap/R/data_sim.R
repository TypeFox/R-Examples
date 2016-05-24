##' @title Simulated binomial data-set over the unit square
##' @description This binomial data-set was simulated by generating a zero-mean Gaussian process over a 30 by 30 grid covering the unit square. The parameters used in the simulation are \code{sigma2=1}, \code{phi=0.15} and \code{kappa=2}. The nugget effect was not included, hence \code{tau2=0}.
##' The variables are as follows:
##' 
##' \itemize{
##'   \item y binomial observations.
##'   \item units.m binomial denominators. 
##'   \item x1 horizontal coordinates.
##'   \item x2 vertical coordinates.
##'   \item S simulated values of the Gaussian process.
##' }
##' 
##' @docType data
##' @keywords datasets
##' @name data_sim
##' @usage data(data_sim)
##' @format A data frame with 900 rows and 5 variables
NULL