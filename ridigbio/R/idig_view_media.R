##' View individual media records.
##' 
##' View all information about a specific media record.
##' @title view media endpoint
##' @param uuid uuid of media record
##' @param ... additional parameters
##' @return nested list of data
##' @author Matthew Collins
##' @export
##'
idig_view_media <- function(uuid, ...){
  view_results <- idig_GET(paste0("view/media/", uuid), ...)
  fmt_view_txt_to_list(view_results)
}

fmt_view_txt_to_list <- function(txt){
  httr::content(txt)
}
