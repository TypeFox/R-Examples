#' A unified framework for input-output operations in R
#'
#' \code{io} provides \code{\link{qread}} for reading in data of various types
#' and \code{\link{qwrite}} for writing data to files of various types.
#' Input or output file types can be inferred from filename extensions or
#' specified explicity.
#' 
#' Use \code{link{io_supported}} to check wehather a data or file type is 
#' supported.
#'
#' Both \code{\link{qread}} and \code{\link{qwrite}} can be readily extended
#' to support additional types by defining specific S3 methods.
#'
#' Additionally, \code{\link{qdraw}} offers a unified interface for plotting
#' to screen or various file formats.
#'
#' @docType package
#' @name io
#' @importFrom grDevices bmp cairo_pdf cairo_ps dev.new dev.off dev.print
#'             jpeg pdf png postscript svg tiff
#' @importFrom utils methods read.table write.table
NULL
