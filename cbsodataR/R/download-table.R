#' Download a table from statistics Netherlands
#' 
#' @param id Identifier of CBS table (can be retrieved from \code{\link{get_table_list}})
#' @param dir Directory where table should be downloaded
#' @param ... Parameters passed on to \code{\link{download_data}}
#' @param cache If metadata is cached use that, otherwise download meta data
#' @param base_url optionally specify a different server. Useful for
#' third party data services implementing the same protocal.
#' 
#' \code{download_table} retrieves all raw meta data and data and stores these as csv
#' files in the directory specified by \code{dir}. It is possible to add a filter. 
#' A filter is specified with \code{<column_name> = <values>} in which \code{<values>} is a character vector.
#' Rows with values that are not part of the character vector are not returned.
#' @export
#' @examples 
#' \dontrun{
#' 
#' # download meta data and data from inflation/Consumer Price Indices
#'  download_table(id="7196ENG")
#' }
download_table <- function(id, ..., dir=id, cache=FALSE, base_url = CBSOPENDATA){
  #TODO add untyped vs typed download
  meta <- download_meta(id=id, dir=dir, cache=cache, base_url = base_url)

  download_data(id=id, path=file.path(dir, "data.csv"), ..., base_url = base_url)
  meta$directory <- dir
  # maybe we should generate a yaml or datapackage.json file?
  invisible(meta)
}

#' @importFrom yaml as.yaml
write_yaml <- function(x, name){
  file_name <- paste0(name, ".yml")
  message("Writing ", file_name, "...")
  writeLines(yaml::as.yaml(x), con=file_name)
}
