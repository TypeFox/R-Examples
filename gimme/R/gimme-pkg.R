#' @name gimme-package
#' @aliases gimme-package
#' @title Group iterative multiple model estimation
#' @description This package contains functions to identify group- and individual-level structural equation models.
#' @author {Stephanie Lane [aut, cre, trl],\cr 
#'    Kathleen Gates [aut], \cr
#'    Peter Molenaar [aut], \cr
#'    Michael Hallquist [ctb], \cr
#'    Hallie Pike [ctb]  
#' }
#' Maintainer: Stephanie Lane \email{slane@@unc.edu}
#' @import lavaan parallel foreach doSNOW doParallel gWidgets2
#' @importFrom grDevices dev.off pdf
#' @importFrom stats aggregate ave complete.cases qchisq qnorm reshape time convolve cov dgamma fft runif sd
#' @importFrom utils capture.output head read.table write.csv write.table file_test read.csv
#' @importFrom MASS ginv
#' @importFrom snow setDefaultClusterOptions makeSOCKcluster
#' @importFrom qgraph qgraph
#' @importFrom igraph graph.adjacency walktrap.community modularity
#' @keywords gimme
#' @details Researchers across varied domains gather multivariate data for each individual unit of study
#' across multiple occasions of measurement. Generally referred to as time series
#' (or in the social sciences, intensive longitudinal) data, examples include
#' psychophysiological processes such as neuroimaging and heart rate variability,
#' daily diary studies, and observational coding of social interactions among dyads.
#' A primary goal for acquiring these data is to understand temporal processes.
#' The gimme package contains several functions for use with these data.
#' These functions include \code{\link{gimmeSEM}}, which provides both group-
#'  and individual-level results by looking across individuals for patterns of
#'   relations among variables. A function that provides group-level results,
#'   \code{\link{aggSEM}}, is included, as well as a function that provides
#'    individual-level results, \code{\link{indSEM}}. The major functions within the gimme package all require the
#'    user to specify the directory containing the data and a directory for output to be stored.
NULL
