#' @title Breast cancer microarray dataset
#'
#' @description 1000 randomly selected probesets from a breast cancer microarray
#' dataset (Farmer et al., 2005).
#'
#' @details These data are the results from an analysis comparing the basal and
#' luminal samples. The apocrine samples are excluded.
#'
#' @format A data frame with 1000 rows and 7 variables:
#' \describe{
#'   \item{ProbeSet:}{Affymetrix probesets from the U133A chip.}
#'   \item{GeneID:}{Gene symbol.}
#'   \item{logFC:}{Log2 fold change for the basal versus luminal comparison.}
#'   \item{AveExpr:}{Mean expression across all samples.}
#'   \item{P.Value:}{P-value for basal versus luminal comparison.}
#'   \item{SD:}{Standard deviation across all samples.}
#'   \item{PCNA.cor:}{Correlation with PCNA (a marker of proliferating cells).}
#' }
#' @references Farmer P, Bonnefoi H, Becette V, Tubiana-Hulin M, Fumoleau P,
#' Larsimont D, Macgrogan G, Bergh J, Cameron D, Goldstein D, Duss S, Nicoulaz
#' AL, Brisken C, Fiche M, Delorenzi M, Iggo R. Identification of molecular
#' apocrine breast tumours by microarray analysis. Oncogene. 2005
#' 24(29):4660-4671.
#'
#' @docType data
#' @name farmer2005
NULL
