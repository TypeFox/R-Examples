#' Example experimental cross data
#'
#' Simulated experimental cross data with some sample mix-ups. The
#' only phenotype is an individual ID. There are 100 individuals
#' genotyped at 1000 markers on 19 autosomes.
#'
#' @docType data
#'
#' @usage data(f2cross)
#'
#' @format An object of class \code{"cross"}. See
#' \code{\link[qtl]{read.cross}} in the R/qtl package for details.
#'
#' @keywords datasets
#'
#' @seealso \code{\link{expr1}}, \code{\link{expr2}}, \code{\link{genepos}}, \code{\link{pmap}}
#'
#' @examples
#' library(qtl)
#' data(f2cross)
#' summary(f2cross)
"f2cross"
