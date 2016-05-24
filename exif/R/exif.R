#' @title Read EXIF data into R
#' @name exif
#' @description exif is a package for reading EXIF media metadata into R,
#' returning it as a list in a similar fashion to jsonlite. It depends on
#' the libexif C library, which must be installed for the package to work.
#' @useDynLib exif
#' @importFrom Rcpp sourceCpp
#' @seealso \code{\link{read_exif}}
#' @docType package
#' @aliases exif exif-package
NULL

#'@title Read EXIF Metadata
#'@description \code{read_exif} reads EXIF metadata from JPEG files,
#'returning it as a data.frame.
#'
#'@param files a vector of files to read in.
#'
#'@return a data.frame, with each row consisting of the metadata for one file in \code{files}. Absent values are
#'represented by an empty string for character columns, and 0 for numeric columns.
#'
#'@examples
#'# A simple example using included images
#'file <- system.file("extdata/dog_test_img.jpg", package="exif")
#'file_metadata <- read_exif(file)
#'
#'@export
read_exif <- function(files){
  return(read_exif_(normalizePath(files, winslash = "\\", mustWork = FALSE)))
}