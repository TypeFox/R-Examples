#' slfm: the sparse latent factor model package for R.
#'
#' Set of tools to find coherent patterns in microarray data
#' using a Bayesian sparse latent factor model (Duarte and Mayrink 2015 -
#' http://link.springer.com/chapter/10.1007%2F978-3-319-12454-4_15).
#' Considerable effort has been put into making slfm fast and memory efficient,
#' turning it an interesting alternative to simpler methods in terms
#' of execution time. It implements versions of the SLFM using both type
#' of mixtures: using a degenerate distribution or a very concentrated
#' normal distribution for the spike part of the mixture. It also implements
#' additional functions to help pre-process the data and fit the model
#' for a large number of arrays. Includes functions to:
#' 
#'   * pre-process a set of matrices
#'   * fit models to a set of matrices
#'   * detailed summary of model fit
#' 
#' @docType package
#' @name slfm-package
NULL