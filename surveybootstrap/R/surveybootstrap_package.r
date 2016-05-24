##' Survey bootstrap variance estimators
##'
##' \code{surveybootstrap} has methods for analyzing data that were collected
##' using network reporting techniques. It includes estimators appropriate for
##' the simple boostrap and the rescaled bootstrap.
##'
##' @docType package
##' @name surveybootstrap
##' @aliases surveybootstrap package-surveybootstrap
##' @import dplyr functional
NULL

##' @useDynLib surveybootstrap
##' @importFrom Rcpp sourceCpp
NULL

##' @importFrom plyr llply
##' @importFrom plyr .
NULL

##' @importFrom stringr str_match
##' @importFrom stringr str_locate
##' @importFrom stringr str_split
NULL

##' @importFrom stats rmultinom
##' @importFrom stats setNames
##' @importFrom stats terms
##' @importFrom stats update
##' @importFrom stats update.formula
##' @importFrom stats xtabs
NULL

##' MU284 population
##'
##' Data used in unit tests for variance estimation.
##' See TODO-Sarndal TODO-sampling package
##' TODO-doc describing unit tests
##'
##' @name MU284
##' @docType data
NULL

