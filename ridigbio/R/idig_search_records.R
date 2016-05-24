##' Function to query the iDigBio API for specimen records
##'
##' Wraps \code{\link{idig_search}} to provide defaults specific to searching 
##' specimen records. Using this function instead of \code{\link{idig_search}} 
##' directly is recommened. 
##' 
##' Queries need to be specified as a nested list structure that will serialize
##' to an iDigBio query object's JSON as expected by the iDigBio API: 
##' \url{https://github.com/iDigBio/idigbio-search-api/wiki/Query-Format}
##' 
##' As an example, the first sample query looks like this in JSON in the API 
##' documentation:
##' \preformatted{
##' {
##'   "scientificname": {
##'     "type": "exists"
##'   },
##'   "family": "asteraceae"
##' }
##' }
##' 
##' To rewrite this in R for use as the rq parameter to 
##' \code{idig_search_records} or \code{idig_search_media}, it would look like
##' this:
##' \preformatted{
##' rq <- list("scientificname"=list("type"="exists"), "family"="asteraceae")
##' }
##' See the Examples section below for more samples of simpler and more complex
##' queries. Please refer to the API documentation for the full functionality
##' availible in queries.
##' 
##' All matching results are returned up to the max_items cap (default 100,000).
##' If more results are wanted, a higher max_items can be passed as an option.
##' This API loads records 5,000 at a time using HTTP so performance with large
##' sets of data is not very good. Expect result sets over 50,000 records to
##' take tens of minutes. You can use the \code{\link{idig_count_records}} or
##' \code{\link{idig_count_media}} functions to find out how many records a
##' query will return; these are fast. 
##' 
##' The iDigBio API will only return 5,000 records at a time but this function 
##' will automatically page through the results and return them all. Limit 
##' and offset are availible if manual paging of results is needed though the 
##' max_items cap still applies. The item count comes from the results header 
##' not the count of actual records in the limit/offset window.
##'
##' Return is a data.frame containing the requested fields (or the default
##' fields). The columns in the data frame are untyped and no factors are pre-
##' built. Attribution and other metadata is attached to the dataframe in the
##' data.frame's attributes. (I.e. \code{attributes(df)})
##' @title Searching of iDigBio records
##' @param rq iDigBio record query in nested list format
##' @param fields vector of fields that will be contained in the data.frame, 
##' limited set returned by default, use "all" to get all indexed fields
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
##' # Simple example of retriving records in a genus:
##' idig_search_records(rq=list(genus="acer"), limit=10)
##' 
##' # This complex query shows that booleans passed to the API are represented
##' # as strings in R, fields used in the query don't have to be returned, and
##' # the syntax for accessing raw data fields:
##' idig_search_records(rq=list("hasImage"="true", genus="acer"), 
##'             fields=c("uuid", "data.dwc:verbatimLatitude"), limit=100)
##' 
##' # Searching inside a raw data field for a string, note that raw data fields
##' # are searched as full text, indexed fields are search with exact matches:
##' 
##' idig_search_records(rq=list("data.dwc:dynamicProperties"="parasite"), 
##'             fields=c("uuid", "data.dwc:dynamicProperties"), limit=100)
##' 
##' # Retriving a data.frame for use with MaxEnt. Notice geopoint is expanded
##' # to two columns in the data.frame: gepoint.lat and geopoint.lon:
##' df <- idig_search_records(rq=list(genus="acer", geopoint=list(type="exists")), 
##'           fields=c("uuid", "geopoint"), limit=10)
##' write.csv(df[c("uuid", "geopoint.lon", "geopoint.lat")], 
##'           file="acer_occurrences.csv", row.names=FALSE)
##'           
##' }
##' @export
##'

idig_search_records <- function(rq, fields=FALSE, max_items=100000, limit=0,
                        offset=0, sort=FALSE, ...) {

  DEFAULT_FIELDS = c('uuid',
                     'occurrenceid',
                     'catalognumber',
                     'family',
                     'genus',
                     'scientificname',
                     'country',
                     'stateprovince',
                     'geopoint',
                     'datecollected',
                     'collector')

  # Validate inputs
  if (!(inherits(rq, "list"))) { stop("rq is not a list") }

  if (!(length(rq) > 0)) { stop("rq must not be 0 length") }

  if (inherits(fields, "logical") && fields == FALSE) {
    fields <- DEFAULT_FIELDS
  }

  if (!(fields == "all" ) && !(inherits(fields, "character"))) {
    stop("Invalid value for fields")
  }

  idig_search(rq=rq, fields=fields, max_items=max_items, limit=limit,
              offset=offset, sort=sort, ...)

}
