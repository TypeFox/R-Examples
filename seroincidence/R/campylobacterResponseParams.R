#' @docType data
#'
#' @name campylobacterResponseParams
#'
#' @title
#' Campylobacter Response Parameters Data
#'
#' @description
#' Monte Carlo sample of longitudinal response parameters \code{A} (peak level)
#' and \code{k} (decay rate) per antibody for Campylobacter collected with SSI ELISA procedure.
#'
#' @usage
#' campylobacterResponseParams
#'
#' @format
#' A list of two dataframes:
#' \describe{
#' \item{\code{A}}{a dataframe containing 4000 peak levels (SSI ELISA units/ml). Named
#' columns contain estimated peak levels for IgG, IgM and IgA antibodies.}
#' \item{\code{k}}{a dataframe containing 4000 decay rates (1/days). Named columns
#' contain estimated decay rates for IgG, IgM and IgA antibodies.}
#' }
#'
#' @references
#' Strid, M. A., Engberg, J., Larsen, L. B., Begtrup, K., Molbak, K., Krogfelt, K. A.\cr
#' "Antibody responses to Campylobacter infections determined by an enzyme-linked immunosorbent assay: 2-year follow-up study of 210 patients"\cr
#' Clinical and Diagnostic Laboratory Immunology 8, no. 2 (March 1, 2001): 314--19. doi:10.1128/CDLI.8.2.314-319.2001.
#'
#' @examples
#'
#' # show first rows of every dataframe contained in campylobacterResponseParams
#' lapply(campylobacterResponseParams, head)
#'
NULL
