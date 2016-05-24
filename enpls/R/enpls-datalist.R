#' Methylalkanes Retention Index Dataset
#'
#' Methylalkanes retention index dataset from Liang et, al.
#'
#' This dataset contains 207 methylalkanes' chromatographic retention index (y) 
#' which have been modeled by 21 molecular descriptors (x).
#'
#' Molecular descriptor types:
#' \itemize{
#' \item Chi path, cluster and path/cluster indices
#' \item Kappa shape indices
#' \item E-state indices
#' \item Molecular electricity distance vector index
#' }
#'
#' @docType data
#' @name alkanes
#' @usage data(alkanes)
#'
#' @format
#' A list with 2 components:
#' \itemize{
#' \item x - data frame with 207 rows (samples) and 21 columns (predictors)
#' \item y - numeric vector of length 207 (response)
#' }
#'
#' @references
#' Yizeng Liang, Dalin Yuan, Qingsong Xu, and Olav Martin Kvalheim. 
#' "Modeling based on subspace orthogonal projections for 
#' QSAR and QSPR research." 
#' Journal of Chemometrics 22, no. 1 (2008): 23--35.
#'
#' @examples
#' data(alkanes)
#' str(alkanes)
NULL
