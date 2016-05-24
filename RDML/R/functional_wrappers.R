#' \code{RDML$AsTable()} wrapper
#' 
#' Read more at \link{RDML.AsTable}
#'   
#' @param obj \code{RDML} object.
#' @param ... \code{AsTable} params.
#' 
#' @name AsTable
#' @rdname astable-function
#' @include RDML.R
#' @export
AsTable <- function(obj, ...) {
  assert_that(is.type(obj, RDML))
  obj$AsTable(...)
}

#' \code{RDML$SetFData()} wrapper
#' 
#' Read more at \link{RDML.SetFData}
#'   
#' @param obj \code{RDML} object.
#' @param ... \code{SetFData} params.
#' 
#' @name SetFData
#' @rdname setfdata-function
#' @include RDML.R
#' @export
SetFData <- function(obj, ...) {
  assert_that(is.type(obj, RDML))
  obj$SetFData(...)
}

#' \code{RDML$GetFData()} wrapper
#' 
#' Read more at \link{RDML.GetFData}
#'   
#' @param obj \code{RDML} object.
#' @param ... \code{GetFData} params.
#' 
#' @name GetFData
#' @rdname GetFData-function
#' @include RDML.R
#' @export
GetFData <- function(obj, ...) {
  assert_that(is.type(obj, RDML))
  obj$GetFData(...)
}

#' \code{RDML$AsDendrogram()} wrapper
#' 
#' Read more at \link{RDML.AsDendrogram}
#'   
#' @param obj \code{RDML} object.
#' @param ... \code{AsDendrogram} params.
#' 
#' @name AsDendrogram
#' @rdname AsDendrogram-function
#' @include RDML.R
#' @export
AsDendrogram <- function(obj, ...) {
  assert_that(is.type(obj, RDML))
  obj$AsDendrogram(...)
}
