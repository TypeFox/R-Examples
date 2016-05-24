#' @title Data to compare the use of blocking and covariate adjustment
#'
#' @description Simulated data to illustrate the effects of blocking
#' versus adjusting for a covariate.
#'
#' @details The experimental manipulation is a new diet versus a
#' standard control diet and the outcome is the amount of food eaten
#' on each diet. Since rats that weigh more at the beginning of the
#' experiment are expected to eat more food, regardless of the diet,
#' it would be beneficial to account for this source of
#' variation. This can be done either through use of blocking or
#' covariate adjustment and data for both designs are included. Note
#' that only one design could be used in a real experiment but here we
#' generate outcome values for two experiments using the same baseline
#' body weight values.
#'
#' For the randomised block design the eight rats are ranked according to
#' baseline body weight and grouped into four blocks of two (the two
#' lightest rats form the first block, the next two the second, and so
#' on). Assignment to treatment group is done within blocks.
#'
#' @format A data frame with 8 rows and 6 variables:
#' \describe{
#' 
#'   \item{weight:}{A baseline measurement of body weight for eight
#' rats. The rows of the data frame are sorted according to weight
#' (lightest to heaviest).}
#' 
#'   \item{block:}{Rats are grouped into four blocks based on body
#' weight.}
#' 
#'   \item{RBD:}{Treatment group when using a randomised block design
#' (RBD).}
#' 
#'   \item{CRD:}{Treatment group when using a completely randomised
#' design (CRD).}
#' 
#'   \item{y.RBD:}{Outcome variable under the RBD.}
#' 
#'   \item{y.CRD:}{Outcome variable under the CRD.}
#' }
#'
#' @examples
#' # Randomised block design
#' summary(aov(y.RBD ~ factor(block) + RBD, data=block.covars))
#'
#' # Completely randomised design with weight as a covariate
#' summary(aov(y.CRD ~ weight + CRD, data=block.covars))
#' 
#' @docType data
#' @name block.covars
NULL
