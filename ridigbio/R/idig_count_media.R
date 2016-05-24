##' Count media records matching a query.
##' 
##' Quickly return a count of the media records matching the query(s) provided.
##' @title Count media endpoint
##' @param rq iDigBio record query in nested list format
##' @param mq iDigBio media query in nested list format
##' @param ... additional parameters
##' @return count of media records matching the query(s)
##' @author Matthew Collins
##' @export
idig_count_media <- function(rq=FALSE, mq=FALSE, ...){
  # This passes through an empty list to get around idig_POST's requirement that
  # rq be present For full API compatability, the post should be completely 
  # empty if the user doesn't specify anything
  query <- list()

  if (inherits(rq, "list") && length(rq) > 0){
    query$rq <- rq
  }

  if (inherits(mq, "list") && length(mq) > 0){
    query$mq <- mq
  }

  view_results <- idig_POST("summary/count/media", body=query, ...)
  fmt_countmedia_txt_to_num(view_results)
}

fmt_countmedia_txt_to_num <- function(txt){
  httr::content(txt)[["itemCount"]]
}
