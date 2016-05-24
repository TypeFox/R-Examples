#' @docType package
#' @description Find Data-only Packages on CRAN
#' @details This package provides a subset of \code{available.packages()}
#' @seealso \code{\link{available_data}}
"_PACKAGE"

#' @title available_data
#' @description Find Data-only Packages on CRAN
#' @details This function returns a data.frame representation of the output of \code{\link[utils]{available.packages}} that contains only known data-only and/or data-heavy packages that are available on CRAN.
#' @param fields A character vector of fields to extract from \code{\link[utils]{available.packages}}.
#' @param \dots Additional arguments passed to \code{\link[utils]{available.packages}}.
#' @return A data.frame containing the requested fields, plus \dQuote{Type} (indicating whether it is a data-only package or a data-supplement to another package) and \dQuote{DataDescription} (providing a short description of the data contained in the package).
#' @examples
#' \dontrun{
#' # get available data packages
#' a <- available_data()
#' }
#' @export
available_data <- 
function(fields = c("Package", "Version", "License"), ...) {
    x <- utils::available.packages(field = fields, ...)
    a <- as.data.frame(x, stringsAsFactors = FALSE)
    kf <- system.file("packages", "packages", package = "crandatapkgs")
    known <- utils::read.csv(kf, stringsAsFactors = FALSE)
    merge(a[, fields], known, all.x = FALSE, all.y = TRUE)
}
