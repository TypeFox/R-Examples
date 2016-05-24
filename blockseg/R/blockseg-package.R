##' blockseg package
##'
##' This package is designed to segment a matrix in blocks with constant values.
##'
##' @section Features: Package for the segmentation of the rows and columns inducing a grid.
##'
##' @section Algorithm: \code{\linkS4class{blockSeg}}, \code{\linkS4class{stab.blockSeg}}
##'
##' @section Technical remarks: Display of the result with \code{\link{plot,blockSeg-method}} and \code{\link{plot,stab.blockSeg-method}} and the evolution
##' with \code{\link{predict,blockSeg-method}} and \code{\link{evolution,stab.blockSeg-method}}.
##'
##' @name blockseg-package
##' @docType package
##' @author Julien Chiquet \email{julien.chiquet@@gmail.com}
##' @author Vincent Brault \email{vincent.brault@@agroparistech.fr}
##'
##' @references Vincent Brault, Julien Chiquet, Celine Levy-Leduc.
##' A Fast Approach for Multiple Change-point Detection in Two-dimensional Data, preprint
##'
##' @importFrom Rcpp sourceCpp
##' @import methods
##' @import Matrix
##' @import ggplot2 
##' @importFrom grDevices gray rgb
##' @importFrom graphics abline axis par title
##' @importFrom utils setTxtProgressBar txtProgressBar
##' @importFrom parallel mclapply detectCores
##' @importFrom reshape2 melt
##' @importFrom stats predict residuals deviance approx rnorm setNames
##' @useDynLib blockseg
NULL
