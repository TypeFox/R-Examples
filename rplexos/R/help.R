#' rplexos: Read and analyze PLEXOS solutions from R
#'
#' Read and analyze PLEXOS solutions from R
#'
#' The following is the typical workflow:
#' 
#' \enumerate{
#'   \item Organize databases: copy all ZIP file solutions into folders with one scenario per folder
#'   \item Process folders: convert the ZIP files with \code{\link{process_folder}}
#'   \item Open databases: create an R object that holds the names of the SQLite databases with \code{\link{plexos_open}}
#'   \item Query data: use the functions documented in \code{\link{query_master}} to extract data
#' }
#' 
#' The list of available properties can be seen with \code{\link{query_property}}.
#'
#' An example of this workflow is available in the "Getting started" vignette.
#' 
#' @docType package
#' @name rplexos
#' @import dplyr
#' @importFrom Rcpp evalCpp
NULL
