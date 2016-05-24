#' sparsebnUtils: Utilities for the sparsebn package.
#'
#' A set of tools for representing and estimating sparse Bayesian
#' networks from continuous and discrete data.
#'
#' This package provides various S3 classes for making it easy to estimate
#' graphical models from data:
#'
#' \itemize{
#' \item \code{\link{sparsebnData}} for managing experimental data with interventions.
#' \item \code{\link{sparsebnFit}} for representing the output of a DAG learning algorithm.
#' \item \code{\link{sparsebnPath}} for representing a solution path of estimates.
#' }
#'
#' The package also provides methods for manipulating these objects and for estimating
#' parameters in graphical models:
#'
#' \itemize{
#' \item \code{\link{estimate.parameters}} for directed graphs.
#' \item \code{\link{estimate.precision}} for undirected graphs.
#' \item \code{\link{estimate.covariance}} for covariance matrices.
#' }
#'
#' Internally, all graph objects may be stored as \code{\link{edgeList}s} (default),
#' or using \code{graphNEL}, \code{igraph}, or \code{network} objects.
#'
#' @docType package
#' @name sparsebnUtils
NULL
