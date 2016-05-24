# Copyright (c) 2016 Nuno Fachada
# Distributed under the MIT License (http://opensource.org/licenses/MIT)

#' micompr: multivariate independent comparison of observations
#'
#' The micompr R package implements a procedure for comparing multivariate
#' samples associated with different groups. The procedure uses principal
#' component analysis to convert multivariate observations into a set of
#' linearly uncorrelated statistical measures, which are then compared using a
#' number of statistical methods. This technique is independent of the
#' distributional properties of samples and automatically selects features that
#' best explain their differences, avoiding manual selection of specific points
#' or summary statistics. The procedure is appropriate for comparing samples of
#' time series, images, spectrometric measures or similar multivariate
#' observations.
#'
#' @references \href{http://arxiv.org/abs/1603.06907}{micompr: An R Package for
#' Multivariate Independent Comparison of Observations}, arXiv preprint.
#'
#' @author Nuno Fachada
#'
#' @note MIT License
#'
#' @seealso \code{\link{cmpoutput}}, \code{\link{micomp}},
#' \code{\link{grpoutputs}}
#'
#' @docType package
#' @name micompr
NULL
#> NULL
