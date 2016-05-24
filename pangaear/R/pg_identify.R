#' Identify information about the Pangaea repository
#'
#' @export
#' @param ... Curl debugging options passed on to \code{\link[httr]{GET}}
#' @return list
#' @examples \dontrun{
#' pg_identify()
#' }
pg_identify <- function(...) {
  res <- oai::id(url = baseoai(), ...)
  structure(as.list(res), class = "pg_identify")
}

#' @export
print.pg_identify <- function(x, ...){
  cat("<Pangaea>", "\n", sep = "")
  cat("  repositoryName: ", x$repositoryName, "\n", sep = "")
  cat("  baseURL: ", x$baseURL, "\n", sep = "")
  cat("  protocolVersion: ", x$protocolVersion, "\n", sep = "")
  cat("  adminEmail: ", x$adminEmail, "\n", sep = "")
  cat("  adminEmail: ", x$adminEmail, "\n", sep = "")
  cat("  earliestDatestamp: ", x$earliestDatestamp, "\n", sep = "")
  cat("  deletedRecord: ", x$deletedRecord, "\n", sep = "")
  cat("  granularity: ", x$granularity, "\n", sep = "")
  cat("  compression: ", paste0(x[ names(x) %in% "compression" ], collapse = ","), "\n", sep = "")
  cat("  description: ", x$description, "\n", sep = "")
}
