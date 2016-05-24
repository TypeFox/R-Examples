#' @include asdreader.r

#' @name soil.asd
#' @title Sample ASD file containing a soil spectrum
#' @docType data
#' @format A binary ASD file
#' @details The spectrum contained in this ASD file was collected
#' on a soil sample from New Zealand using the ASD FieldSpec 3 spectrometer.
#' The file version is ASD 8.0.
#' @examples
#' # Access the location of the ASD file using the following command
#' fn <- asd_file()
#' fn
#' # This function is actually just a shorthand for
#' fn <- system.file("extdata", "soil.asd", package = "asdreader")
#' fn
NULL

#' @name asd_file
#' @title Get location of a sample ASD file
#' @return a character vector storing the location of the sample ASD file
#' @examples
#' fn <- asd_file()
#' fn
asd_file <- function() {
  system.file("extdata", "soil.asd", package = "asdreader")
}
