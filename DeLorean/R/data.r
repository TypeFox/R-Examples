#' Single cell expression data and meta data from Guo et al. (2012).
#' They investigated the expression of 48 genes in 500 mouse embryonic cells.
#'
#' @docType data
#' @keywords datasets
#'
#' @name guo.expr
#' @aliases guo.cell.meta guo.gene.meta
#'
#' @usage data(GuoDeLorean)
#'
#' @format There are three objects in this data:
#' \itemize{
#'   \item guo.expr A matrix of log expression values with
#'     no missing data. Rows are named by genes and columns are
#'     named by cells/samples.
#'   \item guo.gene.meta A data frame containing meta-data
#'     about the genes.
#'   \item guo.cell.meta A data frame containing meta-data
#'     about the cells
#' }
#'
#' @source \url{http://www.sciencedirect.com/science/article/pii/S1534580710001103}
#'
NULL

#' Kouno et al. investigated the transcriptional network controlling how
#' THP-1 human myeloid monocytic leukemia cells differentiate into
#' macrophages. They provide expression values for 45 genes in 960 single
#' cells captured across 8 distinct time points.
#'
#' @docType data
#' @keywords datasets
#'
#' @name kouno.expr
#' @aliases kouno.cell.meta kouno.gene.meta
#'
#' @usage data(KounoDeLorean)
#'
#' @format There are three objects in this data:
#' \itemize{
#'   \item kouno.expr A matrix of log expression values with
#'     no missing data. Rows are named by genes and columns are
#'     named by cells/samples.
#'   \item kouno.gene.meta A data frame containing meta-data
#'     about the genes.
#'   \item kouno.cell.meta A data frame containing meta-data
#'     about the cells
#' }
#'
#' @source \url{http://genomebiology.com/2013/14/10/R118/abstract}
#'
NULL

#' Windram et al. investigated the defense response in Arabidopsis
#' thaliana to the necrotrophic fungal pathogen Botrytis cinerea.
#' They collected data at 24 time points in two conditions for
#' 30336 genes.
#'
#' @docType data
#' @keywords datasets
#'
#' @name windram.expr
#' @aliases windram.cell.meta windram.gene.meta
#'
#' @usage data(WindramDeLorean)
#'
#' @format There are three objects in this data:
#' \itemize{
#'   \item windram.expr A matrix of log expression values with
#'     no missing data. Rows are named by genes and columns are
#'     named by cells/samples.
#'   \item windram.gene.meta A data frame containing meta-data
#'     about the genes.
#'   \item windram.cell.meta A data frame containing meta-data
#'     about the cells
#' }
#'
#' @source \url{http://www.plantcell.org/content/24/9/3530.long}
#'
NULL


## Shalek et al. assayed primary mouse bone-marrow-derived dendritic cells
## under several conditions.
##
## @docType data
## @keywords datasets
##
## @name shalek.A.expr
## @aliases shalek.A.gene.meta shalek.A.cell.meta
##
## @usage data(ShalekDeLorean)
##
## @format There are three objects in this data:
## \itemize{
##   \item shalek.A A matrix of log expression values with
##     no missing data. Rows are named by genes and columns are
##     named by cells/samples.
##   \item shalek.A.gene.meta A data frame containing meta-data
##     about the genes.
##   \item shalek.A.cell.meta A data frame containing meta-data
##     about the cells
## }
##
## @source \url{http://www.nature.com/nature/journal/v510/n7505/full/nature13437.html}
##
NULL
