#'
#' outbreaker: disease outbreak reconstruction using epidemiological and genetic data
#'
#' This package implements the model introduced by Jombart et al. (PLoS Comput. Biol, 2014) for disease outbreak reconstruction using epidemiological and genetic data.
#'
#' Check tutorials and documentation at:
#' \url{https://sites.google.com/site/therepiproject/r-pac/outbreaker}
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @name outbreaker.package
#'
#' @import parallel
#'
#' @importFrom utils packageDescription read.table
#'
#' @importFrom graphics plot barplot lines mtext matplot legend
#'
#' @importFrom grDevices colorRampPalette grey
#'
#' @importFrom graphics arrows boxplot points text
#'
#' @importFrom stats anova density dexp dist lm median rbinom rnorm rpois runif
#'
#' @import igraph
#'
#' @importFrom adegenet num2col transp funky findMutations
#'
#' @importFrom ape dist.dna as.DNAbin
#'
#' @useDynLib outbreaker
#'
NULL
