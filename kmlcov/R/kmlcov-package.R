#' 'artifdata' contains artificial data' obtained with the \code{rnorm} function and should contain 4 cluster of trajectories.
#'
#' @title Artificial data
#' @usage data(artifdata)
#' @name artifdata
#' @docType data
#' @format A [data.frame] containing 500 measures of 50 individuals (trajectories) identified by a
#' column 'id', and the associated taking 'time', 'time2' and 'time3' of some drug for example.
#' In additional there are 2 more colums, 'treatment' is a binary column indicating for example
#' individuals receiving a high dose of some drug and the other receiving a normal dose (coded by 0),
#' 'treatTime' is the 'time' column multiplied by 'treatment'.
#' 
NULL


#'
#' @exportClass Converge
#' @exportClass GlmCluster
#' @exportClass KmlCovList
#' @exportMethod show
#' @exportMethod plot
NULL
