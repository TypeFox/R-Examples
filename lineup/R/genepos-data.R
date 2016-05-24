#' Genomic positions of genes in simulated expression data
#'
#' A table with the genomic positions of genes in the simulated
#' expression data, \code{\link{expr1}} and \code{\link{expr2}}.
#'
#' @docType data
#'
#' @usage data(genepos)
#'
#' @format A data frame with two columns, chromosome and physical position (in Mbp).
#'
#' @keywords datasets
#'
#' @seealso \code{\link{expr1}}, \code{\link{expr2}}, \code{\link{f2cross}}, \code{\link{pmap}}
#'
#' @examples
#' data(genepos)
#'
#' # interplot genetic positions
#' library(qtl)
#' data(pmap)
#' data(f2cross)
#' genepos_interp <- interpPositions(genepos, pmap, pull.map(f2cross))
#' genepos[1:5,] # 'newpos' column is the interpolated cM position
"genepos"
