#' fuzzyforest: an implementation of the fuzzy forest algorithm in R.
#'
#' This package implements fuzzy forests and integrates the fuzzy
#' forests algorithm with the package, \pkg{WGCNA}.
#'
#'
#' @docType package
#' @importFrom randomForest randomForest
#' @importFrom randomForest importance
#' @importFrom randomForest combine
#' @importFrom graphics     plot
#' @importFrom stats        predict
#' @import foreach
#' @import doParallel
#' @import doRNG
#' @import ggplot2
#' @name fuzzyforest
#' @note This work was partially funded by NSF IIS 1251151 and AMFAR 8721SC.
NULL

#' Liver Expression Data from Female Mice
#'
#' A data set containing gene expression levels in liver tissue from female
#' mice. This data set is a subset of the liver expression data set
#' from the WGCNA tutorial \url{http://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/}.
#' The tutorial contains further information about the data set as well as
#' extensive examples of WGCNA.
#'
#' \itemize{
#'   \item The first column contains weight (g) for the 66 mice.
#'   \item The other 3600 columns contain the liver expression levels.
#'  }
#' @docType data
#' @keywords datasets
#' @name Liver_Expr
#' @usage data(Liver_Expr)
#' @format A data frame with 66 rows and 3601
NULL

#' Cardiotocography Data Set
#'
#' A data set containing measurements of fetal heart rate and uterine
#' contraction from cardiotocograms.  This data set was obtained from
#' the [UCI machine learning repository](https://archive.ics.uci.edu/ml/index.html)
#' For our examples we extract a random sub sample of 100 observations.
#'
#' @docType data
#' @keywords datasets
#' @name ctg
#' @usage data(ctg)
#' @format A data frame with 100 rows and 21.
NULL


#' Fuzzy Forest Example
#'
#' An example of a fuzzy_forest object derived from fitting fuzzy forests
#' on the ctg data set.  The source code used to produce example_ff can be
#' seen in the vignette "fuzzyforest_introduction".
#' @docType data
#' @keywords R object
#' @name example_ff
#' @format .RData
NULL


