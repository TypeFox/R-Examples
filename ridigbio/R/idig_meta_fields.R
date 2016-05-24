##' List of fields in iDigBio.
##' 
##' Return a list of media or specimen fields that are contained in iDigBio.
##' @title meta fields endpoint
##' @param type string type of fields to return, defaults to "records"
##' @param subset set of fields to return, "indexed", "raw", or unset for all
##' @param ... additional parameters
##' @return list of fields of the requested type
##' @author Matthew Collins
##' @export
idig_meta_fields <- function(type="records", subset=FALSE, ...){
  # This passes through an empty list to get around idig_POST's requirement that
  # rq be present For full API compatability, the post should be completely 
  # empty if the user doesn't specify anything
  query <- list()

  view_results <- idig_GET(paste0("meta/fields/", type), ...)

  if (subset == "indexed"){
    fmt_metafields_txt_to_indexed(view_results)
  }else if (subset == "raw"){
    fmt_metafields_txt_to_raw(view_results)
  }else{
    fmt_metafields_txt_to_list(view_results)
  }

}

fmt_metafields_txt_to_indexed <- function(txt){
  f <- httr::content(txt)
  f[-match("data", names(f))]
}

fmt_metafields_txt_to_raw <- function(txt){
  httr::content(txt)[["data"]]
}

fmt_metafields_txt_to_list <- function(txt){
  httr::content(txt)
}
