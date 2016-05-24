#' Search for records using keywords/terms to control how your query is interpreted.
#'
#' Returns only those records in which the targeted input is found in
#' association with the specified search terms.
#'
#' @export
#'
#' @param specificepithet Taxonomic specific epithet, e.g. (sapiens in Homo sapiens) (character)
#' @param genus Taxonomic genus (character)
#' @param family Taxonomic family (character)
#' @param order Taxonomic order (character)
#' @param class Taxonomic class (character)
#' @param limit (numeric) Limit on the number of records returned. If >1000 results, we use
#' a cursor internally, but you should still get up to the results you asked for. See also 
#' \code{\link{bigsearch}} to get larger result sets in a text file via email.
#' @param compact Return a compact data frame (logical)
#' @param year Year (numeric) or range of years designated by comparison
#'  operators "<", ">", "<=" or ">=". You can pass in more than one of these
#'  queries, in a vector. See example below. (character)
#' @param date Event date associated with this occurrence record; yyyy-mm-dd
#'  or the range yyyy-mm-dd/yyyy-mm-dd (character)
#' @param mappable Record includes valid coordinates in decimal latitude and
#'  decimal longitude (logical)
#' @param error Coordinate uncertainty in meters (numeric) or range of uncertainty
#'  values designated by comparison operators "<", ">", "<=" or ">=" (character)
#' @param continent Continent to search for occurrence (character)
#' @param cntry Country to search for occurrence (character)
#' @param stateprovince State or province to search for occurrence (character)
#' @param county County to search for occurrence (character)
#' @param island Island to search for occurrence (character)
#' @param igroup Island group to search for occurrence (character)
#' @param inst Code name for the provider/institution of record (character)
#' @param id Provider's unique identifier for this occurrence record (character)
#' @param catalognumber Provider's catalog number or other ID for this record (character)
#' @param collector Collector name (character)
#' @param type Type of record; "specimen" or "observation" (character)
#' @param hastypestatus Specimen associated with this record is identified as a
#'    holotype, paratype, neotype, etc. (character)
#' @param media Record also references associated media, such as a film or video (logical)
#' @param rank TBD (numeric)
#' @param tissue Record is likely to reference tissues (logical)
#' @param resource Identifier for the resource/dataset from which the record was
#'  indexed (character)
#' @param verbose Print progress and information messages. Default: TRUE
#' @param ... Curl arguments passed on to \code{\link[httr]{GET}}
#'
#' @return A data frame of search results
#' @references \url{https://github.com/VertNet/webapp/wiki/The-API-search-function}
#' @details \code{\link{searchbyterm}} builds a query from input parameters based on
#'    Darwin Core (dwc) terms (for the full list of terms, see
#'    \url{https://code.google.com/p/darwincore/wiki/DarwinCoreTerms}).
#'    The query string is appended to the base URL for VertNet
#'    search requests. View the query string for specification of
#'    dwc terms used in the search.
#'
#' @examples \dontrun{
#' # Limit the number of records returned to <1000; use bigsearch() for >1000 records
#' (out <- searchbyterm(class = "aves", state = "california", limit = 10))
#'
#' # Find multiple species
#' (out <- searchbyterm(genus = "ochotona", specificepithet = "(princeps OR collaris)", limit=10))
#'
#' # Specifying a single year (no quotes) or range of years (use quotes)
#' (out <- searchbyterm(class = "aves", state = "california", year = 1976, limit = 10))
#' (out <- searchbyterm(class = "aves", state = "california", year = ">=1976", limit = 10))
#'
#' # Specifying a range (in meters) for uncertainty in spatial location (use quotes)
#' out <- searchbyterm(class = "aves", state = "nevada", error = "<25")
#' out <- searchbyterm(class = "aves", state = "california", year = 1976, error = "<=1000")
#'
#' # Specifying records by event date (use quotes)
#' out <- searchbyterm(class = "aves", state = "california", date = "2009-03-25")
#' # ...but specifying a date range may not work
#' out <- searchbyterm(specificepithet = "nigripes", date = "1935-09-01/1935-09-30")
#'
#' # Pass in curl options for curl debugging
#' library("httr")
#' out <- searchbyterm(class = "aves", limit = 10, config=verbose())
#' # out <- searchbyterm(class = "aves", limit = 500, config=timeout(1))
#'
#' # Use more than one year query
#' searchbyterm(class = "aves", year = c(">=1976", "<=1986"))
#' }

searchbyterm <- function(specificepithet = NULL, genus = NULL, family = NULL, order = NULL,
                  class = NULL, limit = 1000, compact = TRUE, year = NULL, date = NULL,
                  mappable = NULL, error = NULL, continent = NULL, cntry = NULL,
                  stateprovince = NULL, county = NULL, island = NULL, igroup = NULL,
                  inst = NULL, id = NULL, catalognumber = NULL, collector = NULL, type = NULL,
                  hastypestatus = NULL, media = NULL, rank = NULL, tissue = NULL,
                  resource = NULL, verbose = TRUE, ...) {

  args <- compact(list(specificepithet = specificepithet, genus = genus, family = family,
        order = order, class = class, eventdate = date,
        mappable = ab(mappable), coordinateuncertaintyinmeters = error,
        continent = continent, country = cntry, stateprovince = sp_paren(stateprovince),
        county = county, island = island, islandgroup = igroup,
        institutioncode = inst, occurrenceid = id, catalognumber = catalognumber,
        recordedby = collector, type = type, hastypestatus = hastypestatus,
        media = ab(media), rank = rank, tissue = ab(tissue), resource = resource))
  args <- compact(c(args, combyr(year)))
  vertwrapper(fxn = "searchbyterm", args = args, lim = limit, compact = compact, verbose = verbose, ...)
}

sp_paren <- function(x) {
  if (!is.null(x)) {
    if (grepl("\\sAND\\s|\\sOR\\s|\\sNOT\\s", x)) {
      sprintf("(%s)", x)
    } else {
      x
    }
  } else {
    NULL
  }
}

ab <- function(x){
  if (is.null(x)) {
    NULL
  } else {
    if (x) 1 else 0
  }
}
