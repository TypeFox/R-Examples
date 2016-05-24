#' Statistical Power Simulation for Testing the Rasch Model
#' 
#' Statistical power simulation for testing the Rasch Model 
#' based on a three-way analysis of variance design with mixed classification.
#' 
#' @docType package
#' @name pwrRasch
#' 
#' @importFrom graphics plot legend points axis box lines mtext
#' @importFrom stats pf printCoefmat rnorm runif aov
#' 
#' @author 
#'  Takuya Yanagida [aut,cre] \email{takuya.yanagida@@univie.ac.at},
#'  Jan Steinfeld [aut] \email{jan.steinfeld@@univie.ac.at},
#'  Thomas Kiefer [ctb]
#'  
#'  Maintainer: Takuya Yanagida <takuya.yanagida@@univie.ac.at>
#'  
#' @seealso \link{aov.rasch}, \link{pwr.rasch}
#' 
#' @references 
#' Kubinger, K. D., Rasch, D., & Yanagida, T. (2009). On designing data-sampling for Rasch model calibrating an achievement test. \emph{Psychology Science Quarterly, 51}, 370-384.
#'
#' Kubinger, K. D., Rasch, D., & Yanagida, T. (2011). A new approach for testing the Rasch model. \emph{Educational Research and Evaluation, 17}, 321-333.
#' 
#' Verhelst, N. D. (2008). An efficient MCMC algorithm to sample binary matrices with fixed marginals. \emph{Psychometrika, 73}(4), 705-728.
#' 
#' Verhelst, N., Hatzinger, R., & Mair, P. (2007). The Rasch sampler. \emph{Journal of Statistical Software, 20}(4), 1-14.
#'
NULL
#' Sample of test data from subtest 2 of the Adaptive Intelligence Diagnosticum (AID3; Kubinger \& Holocher-Ertl, 2014)
#'
#' A dataset containing the test data of 300 childen (drawn randomly from the original dataset).
#' The variables are as follows:
#'
#' @format A data frame with 300 rows and 28 variables:
#' \itemize{
#'   \item ID: ID variable of each testee
#'   \item age_in_month: the age of the testperson in month
#'   \item sex: gender of the testee
#'   \item country: country of the testee
#'   \item stage: stage of the data collection
#'   \item it1...it18: items of the subtest 2
#' }
"aid_st2"
