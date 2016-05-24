##' Function to query the iDigBio API for media records
##' 
##' Also see \code{\link{idig_search_records}} for the full examples of all the 
##' parameters related to searching iDigBio.
##' 
##' Wraps \code{\link{idig_search}} to provide defaults specific to searching 
##' media records. Using this function instead of \code{\link{idig_search}} 
##' directly is recommened. Record queries and media queries objects are allowed
##' (rq and mq parameters) and media records returned will match the 
##' requirements of both.
##' 
##' This function defaults to returning all indexed media record fields. 
##' @title Searching of iDigBio media records
##' @param mq iDigBio media query in nested list format
##' @param rq iDigBio record query in nested list format
##' @param fields vector of fields that will be contained in the data.frame, 
##' defaults to "all" which is all indexed fields
##' @param max_items maximum number of results allowed to be retrieved (fail
##' -safe)
##' @param limit maximum number of results returned
##' @param offset number of results to skip before returning results
##' @param sort vector of fields to use for sorting, UUID is always appended to 
##' make paging safe
##' @param ... additional parameters
##' @return a data frame
##' @author Matthew Collins
##' @examples
##' \dontrun{
##' # Searching for media using a query on related specimen information - first 
##' # 10 media records with image URIs related to a specimen in the genus Acer:
##' df <- idig_search_media(rq=list(genus="acer"), 
##'                         mq=list("data.ac:accessURI"=list("type"="exists")), 
##'                         fields=c("uuid","data.ac:accessURI"), limit=10)
##' }
##' @export
##'
idig_search_media <- function(mq=FALSE, rq=FALSE, fields=FALSE, 
                              max_items=100000, limit=0, offset=0, sort=FALSE,
                              ...) {

  DEFAULT_FIELDS = "all"

  # Validate inputs
  #if (!(inherits(rq, "list"))) { stop("rq is not a list") }

  if (!(length(rq) > 0) && !(length(rq) > 0)) {
    stop("mq or rq must not be 0 length")
  }

  if (inherits(fields, "logical") && fields == FALSE) {
    fields <- DEFAULT_FIELDS
  }

  if (!(fields == "all" ) && !(inherits(fields, "character"))) {
    stop("Invalid value for fields")
  }

  idig_search(type="media", mq=mq, rq=rq, fields=fields, max_items=max_items, 
              limit=limit, offset=offset, sort=sort, ...)

}
