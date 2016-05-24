#' @title Valproic acid (VPA) data set
#'
#' @description A mouse split-plot experiment with two treatments.
#'
#' @details Six pregnant female mice were randomly assigned to receive
#' an injection of valproic acid (n = 3) or saline (n = 3). The
#' offspring of these mice (n = 24) were then randomly assigned to
#' receive an injection of the glutamate receptor antagonist MPEP (n =
#' 12) or saline (n = 12). There are two levels of randomisation: the
#' pregnant females and their offspring.
#'
#' This is a subset of the full data set from Mehta et al., which
#' contained fourteen litters. The litters were selected to have an
#' equal sample size across all conditions for illustrative
#' purposes. The complete data can be found in the supplementary
#' material of Lazic and Essioux (2013).
#'
#' @format A data frame with 24 rows and 5 variables:
#' \describe{
#'   \item{litter:}{Mice were from one of six litters.}
#'   \item{sex:}{Male or Female.}
#'   \item{group:}{VPA or saline (SAL).}
#'   \item{drug:}{MPEP or saline (SAL).}
#'   \item{activity:}{Locomotor activity.}}
#'
#' @references Mehta MV, Gandal MJ, Siegel SJ (2011). mGluR5-antagonist mediated reversal of elevated stereotyped, repetitive behaviors in the VPA model of autism. \emph{PLoS ONE} 6(10):e26077.
#' 
#' @references Lazic SE, Essioux L (2013). Improving basic and
#' translational science by accounting for litter-to-litter variation
#' in animal models. \emph{BMC Neuroscience} 14:37.
#'
#' 
#' @docType data
#' @name VPA
NULL
