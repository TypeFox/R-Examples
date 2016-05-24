#' UPB data
#'
#' Data from a survey study that was part of the Interdisciplinary Project for the Optimization of Separation trajectories (IPOS). This large-scale project involved the recruitment of individuals who divorced between March 2008 and March 2009 in four major courts in Flanders. It aimed to improve the quality of life in families during and after the divorce by translating research findings into practical guidelines for separation specialists and by promoting evidence-based policy.
#' This dataset involves a subsample of 385 individuals, namely those who responded to a battery of questionnaires related to romantic relationship and breakup characteristics (De Smet, 2012).
#' @format A data frame with 385 rows and 9 variables:
#' \describe{
#'   \item{att}{self-reported anxious attachment level (standardized)}
#'   \item{attbin}{binary version of self-reported anxious attachment level: 1 = higher than sample mean, 0 = lower than sample mean}
#'   \item{attcat}{multicategorical version of self-reported anxious attachment level: \code{L} = low, \code{M} = intermediate, \code{H} = high}
#'   \item{negaff}{level of self-reported experienced negative affectivity (standardized)}
#   \item{negaffbin}{binary version of self-reported experienced negative affectivity: 1 = higher than sample mean, 0 = lower than sample mean}
#   \item{negaffcat}{multicategorical version of level of self-reported experienced negative affectivity: \code{L} = low, \code{M} = intermediate, \code{H} = high}
#'   \item{initiator}{initiator of the divorce}
#'   \item{gender}{gender: \code{F} = female, \code{M} = male}
#'   \item{educ}{education level: either \code{H} = high (at least a bachelor's degree), \code{M} = intermediate (having finished secondary school) or \code{L} = low (otherwise)}
#'   \item{age}{age (in years)}
#'   \item{UPB}{binary variable indicating whether the individual reported having displayed unwanted pursuit behavior(s) towards the ex-partner}
#' }
#' @source
#' Ghent University and Catholic University of Louvain (2010). \emph{Interdisciplinary Project for the Optimisation of Separation trajectories - divorce and separation in Flanders}. 
#' \cr\url{http://www.scheidingsonderzoek.ugent.be/index-eng.html}
#' @references
#' De Smet, O., Loeys, T., & Buysse, A. (2012). Post-Breakup Unwanted Pursuit: A Refined Analysis of the Role of Romantic Relationship Characteristics. \emph{Journal of Family Violence}, \bold{27}(5), 437-452.
#' 
#' Loeys, T., Moerkerke, B., De Smet, O., Buysse, A., Steen, J., & Vansteelandt, S. (2013). Flexible Mediation Analysis in the Presence of Nonlinear Relations: Beyond the Mediation Formula. \emph{Multivariate Behavioral Research}, \bold{48}(6), 871-894.
#' @name UPBdata
NULL