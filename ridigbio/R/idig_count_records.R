##' Count specimen records matching a query.
##' 
##' Quickly return a count of the specimen records matching the query(s) 
##' provided.
##' @title Count record endpoint
##' @param rq iDigBio record query in nested list format
##' @param ... additional parameters
##' @return count of specimen records matching the query(s)
##' @author Matthew Collins
##' @export
idig_count_records <- function(rq=FALSE, ...){
  # This passes through an empty list to get around idig_POST's requirement that
  # rq be present For full API compatability, the post should be completely 
  # empty if the user doesn't specify anything
  query <- list()

  if (inherits(rq, "list") && length(rq) > 0){
    query$rq <- rq
  }

  view_results <- idig_POST("summary/count/records", body=query, ...)
  fmt_countrecords_txt_to_num(view_results)
}

fmt_countrecords_txt_to_num <- function(txt){
  httr::content(txt)[["itemCount"]]
}
