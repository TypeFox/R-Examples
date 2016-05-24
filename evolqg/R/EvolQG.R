#' EvolQG
#'
#' @name evolqg
#' @docType package
#' @import plyr
#' @importFrom graphics abline arrows axis layout mtext par plot text
#' @importFrom grDevices rgb
#' @importFrom methods is
#' @importFrom stats anova confint cor cor.test cov cov2cor lm mahalanobis princomp quantile reorder residuals rnorm runif sd var
#' @importFrom utils write.csv write.table
NULL

#' Example multivariate data set
#'
#' Simulated example of 4 continuous bone lengths from 5 species.
#'
#' \itemize{
#' \item humerus 
#' \item ulna 
#' \item femur 
#' \item tibia 
#' \item species 
#' }
#'
#' @docType data
#' @keywords datasets
#' @name dentus
#' @usage data(dentus)
#' @format A data frame with 300 rows and 5 variables
NULL

#' Tree for dentus example species
#'
#' Hypothetical tree for the species in the dentus data set.
#'
#' @docType data
#' @keywords datasets
#' @name dentus.tree
#' @usage data(dentus.tree)
#' @format ape tree object
NULL