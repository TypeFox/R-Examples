#' Deprecated functions in package rbefdata
#'
#' These following functions are renamed or removed. Please use new functions instead.
#'
#' \tabular{ll}{
#'  \strong{Deprecated functions} \tab \strong{Replacement} \cr
#'  \code{bef.portal.get.dataset_list} \tab \code{\link{bef.portal.get.datasets.for_keyword}} \cr
#'  \code{bef.portal.get.proposal} \tab \code{\link{bef.portal.get.datasets.for_proposal}} \cr
#' }
#'
#' @name rbefdata-deprecated
#' @aliases bef.portal.get.dataset_list bef.portal.get.proposal
#' @export bef.portal.get.dataset_list bef.portal.get.proposal
#' @keywords internal

## version 0.3 -> 0.4 deprecates
#' @rdname rbefdata-deprecated
bef.portal.get.dataset_list <- function(...) {
  .Deprecated("bef.portal.get.datasets_for_keyword")
  bef.portal.get.datasets.for_keyword(...)
}

#' @rdname rbefdata-deprecated
bef.portal.get.proposal <- function(...) {
  .Deprecated("bef.portal.get.datasets_for_proposal")
  bef.portal.get.datasets.for_proposal(...)
}

