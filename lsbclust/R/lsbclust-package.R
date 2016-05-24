#' @name lsbclust-package
#' @aliases lsbclust-package
#' @title Least Squares Latent Class Matrix Factorization
#' @docType package
#' @description Funtions for least squares latent class matrix factorizations.
#' @author Pieter C. Schoonees [aut, cre], Patrick J.F. Groenen [aut]
#' @references Van Rosmalen, J., Van Herk, H., & Groenen, P. J. F. (2010). Identifying response styles: A latent-class bilinear multinomial logit model. \emph{Journal of Marketing Research}, 47(1), 157-172.
#' @keywords package
#' @import ggplot2
#' @importFrom plyr alply mapvalues
#' @importFrom clue solve_LSAP cl_class_ids is.cl_partition is.cl_hard_partition as.cl_class_ids cl_agreement
#' @importFrom gridExtra grid.arrange
#'@importFrom Rcpp evalCpp
#'@useDynLib lsbclust
NULL