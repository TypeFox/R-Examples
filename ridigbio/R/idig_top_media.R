##' Top media records summaries.
##' 
##' Summarize the count of media records in iDigBio according to unique values
##' in the fields passed. This operates similarly to a SELECT DISTINCT count(
##' field_name) query in SQL. When multiple fields are passed, the summaries are
##' nested eg fields=c("country", "genus") would result in counting the top 10
##' genera in each of the top 10 countries for a total of 100 counts.
##' @title Top media endpoint
##' @param rq iDigBio record query in nested list format
##' @param mq iDigBio media query in nested list format
##' @param top_fields vector of field names to summarize by
##' @param count maximum number of results to return, capped at 1000
##' @param ... additional parameters
##' @return nested list of field values with counts of media records
##' @author Matthew Collins
##' @export
##' 
idig_top_media <- function(rq=FALSE, mq=FALSE, top_fields=FALSE, count=0, ...){
  # This passes through an empty list to get around idig_POST's requirement that
  # rq be present. For full API compatability, the post should be completely 
  # empty if the user doesn't specify anything
  query <- list()

  if (inherits(rq, "list") && length(rq) > 0){
    query$rq <- rq
  }

  if (inherits(mq, "list") && length(mq) > 0){
    query$mq <- mq
  }

  if (inherits(top_fields, "character") && length(top_fields) > 0){
    query$top_fields <- top_fields
  }

  if (count > 0){
    query$count <- count
  }

  view_results <- idig_POST("summary/top/media", body=query, ...)
  fmt_topmedia_txt_to_list(view_results)
}

fmt_topmedia_txt_to_list <- function(txt){
  httr::content(txt)
}
