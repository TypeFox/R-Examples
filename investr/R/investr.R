##' investr: a package for inverse estimation in R
##' 
##' Inverse estimation, also referred to as the calibration problem, is a 
##' classical and well-known problem in regression. In simple terms, it involves 
##' the use of an observed value of the response (or specified value of the mean 
##' response) to make inference on the corresponding unknown value of the 
##' explanatory variable. 
##'
##' A detailed introduction to investr has been published in The R Journal: 
##' "investr: An R Package for Inverse Estimation", 
##' \url{http://journal.r-project.org/archive/2014-1/greenwell-kabban.pdf}. You 
##' can track development at \url{https://github.com/w108bmg/investr}. To report 
##' bugs or issues, contact the main author directly or submit them to 
##' \url{https://github.com/w108bmg/investr/issues}. 
##'
##' As of right now, \code{investr} supports (univariate) inverse estimation 
##' with objects of class:
##' \itemize{
##'   \item{\code{lm}} --- linear models (multiple predictor variables allowed)
##'   \item{\code{glm}} --- generalized linear models (multiple predictor variables allowed)
##'   \item{\code{nls}} --- nonlinear least-squares models
##'   \item{\code{lme}} --- linear mixed-effects models (fit using the 
##'     \code{nlme} package)
##' }
##' 
##' @docType package
##' @name investr
NULL
