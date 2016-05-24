##' Base function to query the iDigBio API
##'
##' This function is wrapped for media and specimen record searches. Please
##' consider using \code{\link{idig_search_media}} or
##' \code{\link{idig_search_records}} instead as they supply nice defaults to
##' this function depending on the type of records desired.
##'
##' Fuller documentation of parameters is in the
##' \code{\link{idig_search_records}} function's help.
##'
##' Exported to facilitate wrapping this package in other packages.
##' @title Basic searching of iDigBio records
##' @param type string type of records to query, defaults to "records"
##' @param mq iDigBio media query in nested list format
##' @param rq iDigBio record query in nested list format
##' @param fields vector of fields that will be contained in the data.frame
##' @param max_items maximum number of results allowed to be retrieved (fail
##' -safe)
##' @param limit maximum number of results returned
##' @param offset number of results to skip before returning results
##' @param sort vector of fields to use for sorting, UUID is always appended to
##' make paging safe
##' @param ... additional parameters
##' @return a data frame
##' @author Francois Michonneau
##' @export
##' @examples
##' \dontrun{
##' # Ten media records related to genus Acer specimens
##' idig_search(type="media", rq=list(genus="acer"), limit=10)
##' }
##'
idig_search <- function(type="records", mq=FALSE, rq=FALSE, fields=FALSE,
                        max_items=100000, limit=0, offset=0, sort=FALSE, ...) {
  # Construct body of request to API
  query <- list(offset=offset)

  # Force sorting by UUID so that paging will be reliable ie the 25,000th item
  # is always the 25,000th item even when requesting the 6th page. This
  # has been benchmarked and appears make things ~20% slower on a gigabit
  # connection: 66s for 100,000 limit
  if (!inherits(sort, "logical")) {
    query[["sort"]] <- c(sort, "uuid")
  }else{
    query[["sort"]] <- c("uuid")
  }

  if (!inherits(rq, "logical")) {
    query$rq=rq
  }

  if (!inherits(mq, "logical")) {
    query$mq=mq
  }

  # Adjust fields to request from the API
  field_lists <- build_field_lists(fields, type)
  fields <- field_lists$fields
  query <- append(query, field_lists$query)

  if (limit > 0){
    query$limit <- limit
  }else{
    query$limit <- max_items # effectivly use iDigBio's max page size
  }

  # tricks to get inside loop first time
  m <- matrix(nrow=0, ncol=length(fields))
  #    res <- data.frame(res, stringsAsFactors = FALSE)
  dat <- data.frame(m, stringsAsFactors=FALSE)
  colnames(dat) <- fields
  item_count <- 1

  # loop until we either have all results or all results the user wants
  while (nrow(dat) < item_count && (limit == 0 || nrow(dat) < limit)){
    search_results <- idig_POST(paste0("search/", type), body=query, ...)
    #print(paste0(Sys.time(), " completed query"))
    # Slight possibility of the number of items changing as we go due to inserts
    # deletes at iDigBio, put this inside the loop to keep it current
    item_count <- fmt_search_txt_to_itemCount(search_results)
    if ((limit == 0 || limit > max_items) && item_count > max_items){
      stop(paste0("Search would return more than ", max_items,
                  " results. See max_items argument."))
    }

    dat <- plyr::rbind.fill(dat, fmt_search_txt_to_df(search_results, fields))

    query$offset <- nrow(dat)
    if (limit > 0){
      query$limit <- limit - nrow(dat)
    }
  }

  # Metadata as attributes on the df
  a <- attributes(dat)
  a[["itemCount"]] <- item_count
  a[["attribution"]] <- fmt_search_txt_to_attribution(search_results)
  attributes(dat) <- a

  dat
}

fmt_search_txt_to_attribution <- function(txt){
  httr::content(txt)[["attribution"]]
}

fmt_search_txt_to_itemCount <- function(txt){
  httr::content(txt)$itemCount
}

fmt_search_txt_to_df <- function(txt, fields) {
  # Check returned results for common errors
  if (!exists("items", httr::content(txt))){
    stop("Returned results do not contain any content")
  }

  #Before continuing to add error handling, let's settle on a pattern.

  search_items <- jsonlite::fromJSON(httr::content(txt, as="text"))[["items"]]
  res <- data.frame(search_items[["indexTerms"]], 
                    search_items[["data"]], stringsAsFactors=FALSE)

  # Append "data." to the data field names. Also, for some reason ":" gets 
  # changed to "." in the data field names when making the df. Need the if
  # statements because helpfully paste0("str", NULL) => "str".
  n <- c()
  if (length(names(search_items[["indexTerms"]])) > 0){
    n <- c(n, names(search_items[["indexTerms"]]))
  }
  if (length(names(search_items[["data"]])) > 0){
    n <- c(n, paste0("data.", names(search_items[["data"]])))
  }
  colnames(res) <- n
  
  # Fixup geopoint into two fields for convenience
  # Doing this inside here because the rbind.fill function seems to pack
  # list fields into nested lists in the last record of the first df made. It's 
  # weird. Would be nicer to do this outside the paging loop otherwise.
  if ("geopoint" %in% colnames(res)){
    res[["geopoint.lon"]] <- res[["geopoint"]][[1]]
    res[["geopoint.lat"]] <- res[["geopoint"]][[2]]
    res$geopoint <- NULL
  }

  res
}


##' Build fields and fields_exclude for queries.
##'
##' Given the desired fields to be returned, intelligently add an exclusion for
##' the data array if warranted and handle the "all" keyword. And do so without
##' setting both fields and fields_exclude due to fact that the API will return
##' wrong results if are passed. This is still posssible if the user
##' deliberately sets both. Not exported.
##' @param fields character vector of fields user wants returned
##' @param type type of records to get fields for
##' @return list list with fields key for df fields and query key for parameters
##' to be merged with the query sent
build_field_lists <- function(fields, type) {
  ret <- list()
  ret$query = list()
  # Here Alex says to eat "all" rather than pass it through to the API
  if (inherits(fields, "character") && fields != "all" && length(fields) > 0 ){
    ret$fields <- fields
    ret$query$fields <- fields
  } else {
    # When a field parameter is passed then the un-requested raw data is
    # already dropped because it's not a requested field. When no field
    # parameter is passed then drop it manually since by default we will not
    # return data fields and this saves significant transfer.
    ret$query$fields_exclude <- "data"
    # Load up all fields possible
    ret$fields <- names(idig_meta_fields(type=type, subset="indexed"))
  }
  
  # Fixup geopoint into two fields. There is also a parallel fixup inside the
  # fmt_search_txt_to_df() function. Preserve field order that the user 
  # specified.
  if ("geopoint" %in% ret[["fields"]]){
    i <- match("geopoint", ret[["fields"]])
    ret[["fields"]][[i]] <- "geopoint.lon"
    ret[["fields"]] <- append(ret[["fields"]], "geopoint.lat", i)
  }
  ret
}
