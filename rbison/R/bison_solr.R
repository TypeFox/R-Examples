#' Search for and collect occurrence data from the USGS Bison API using their solr endpoint.
#'
#' This fxn is somewhat similar to \code{\link{bison}}, but interacts with the SOLR
#' interface \url{http://bison.usgs.ornl.gov/doc/api.jsp#solr} instead of the OpenSearch interface
#' \url{http://bison.usgs.ornl.gov/doc/api.jsp#opensearch}, which \code{\link[rbison]{bison}} uses.
#'
#' @export
#'
#' @param decimalLatitude Geographic coordinate that specifies the north south position
#' of a location on the Earth surface.
#' @param decimalLongitude	Geographic coordinate that specifies the east-west position
#' of a location on the Earth surface.
#' @param year The year the collection was taken.
#' @param providerID (character) Unique identifier assigned by GBIF.
#' @param resourceID (character) A unique identifier that is a concatentation of the provider
#' identifier and the resource id seperated by a comma.
#' @param pointPath	A dynamic field that contains the location in longitude and
#' latitude followed by the basis of record and an optional Geo (Spatial) precision.
#' Geo (Spatial) precision is an added descriptor when the record is a county centroid.
#' @param basisOfRecord	One of these enumerated values: Observation, Germplasm,
#' Fossil, Specimen, Literature, Unknown, or Living.
#' @param eventDate	The date when the occurrence was recorded.
#' @param computedCountyFips County FIPS code conforming to standard FIPS 6-4 but with leading
#' zeros removed.
#' @param computedStateFips The normalized case sensitive name. For example
#' q=state_code:"New Mexico" will return all of the occurrences from New Mexico.
#' @param scientificName	The species scientific name that is searchable in a case
#' insensitive way.
#' @param hierarchy_homonym_string hierarachy of the accepted or valid species name starting at
#' kingdom. If the name is a taxonomic homonym more than one string is provided seperated by ';'.
#' @param TSNs Accepted or valid name is provided. If the name is a taxonmic homonym more
#' than one TSN is provided.
#' @param recordedBy Individual responsible for the scientific record.
#' @param occurrenceID Non-persistent unique identifier.
#' @param catalogNumber Unique key for every record (occurrence/row) within a dataset that is not
#' manipulated nor changed (nor generated, if not provided) during the data ingest.
#' @param ITIScommonName Common name(s) from ITIS, e.g. "Canada goose"
#' @param kingdom Kingdom name, from GBIF raw occurrence or BISON provider.
#' @param callopts Further args passed on to httr::GET for HTTP debugging/inspecting. In
#' \code{bison}, \code{bison_providers}, and \code{bison_stats}, \code{...} is used instead of
#' callopts, but \code{...} is used here to pass additional Solr params.
#' @param ... Additional SOLR query arguments. See details.
#' @param verbose Print message with url (TRUE, default).
#'
#' @return An object of class bison_solr - which is a list with slots for number of
#' records found (num_found), records, highlight, or facets.
#' @details Some SOLR search parameters:
#' \itemize{
#'  \item{fl} {Fields to return in the query}
#'  \item{rows} {Number of records to return}
#'  \item{sort} {Field to sort by, see examples}
#'  \item{facet} {Facet or not, logical}
#'  \item{facet.field} {Fields to facet by}
#' }
#'
#' You can also use highlighting in solr search, but I'm not sure I see a use case for it
#' with BISON data, though you can do it with this function.
#'
#' For a tutorial see here \url{http://lucene.apache.org/solr/3_6_2/doc-files/tutorial.html}
#' @seealso \code{\link{bison_tax}} \code{\link{bison}}
#'
#' The USGS BISON Solr installation version as of 2014-10-14 was 4.4.
#'
#' @examples \dontrun{
#' bison_solr(scientificName='Ursus americanus')
#' }
#'
#' @examples \dontrun{
#' bison_solr(scientificName='Ursus americanus', computedStateFips='New Mexico',
#'  fl="scientificName")
#'
#' bison_solr(scientificName='Ursus americanus', computedStateFips='New Mexico',
#'  rows=50, fl="eventDate,scientificName")
#'
#' bison_solr(providerID = 220)
#' bison_solr(resourceID = '220,200080')
#'
#' bison_solr(eventDate = '2010-08-08T00:00Z')
#'
#' bison_solr(TSNs = 174773)
#'
#' bison_solr(occurrenceID = 576630651)
#'
#' bison_solr(catalogNumber = 'OBS101299944')
#'
#' bison_solr(ITIScommonName = "Canada goose")
#'
#' bison_solr(kingdom = "Animalia")
#' bison_solr(kingdom = "Plantae")
#'
#' # Mapping
#' out <- bison_solr(scientificName='Ursus americanus', rows=200)
#' bisonmap(out)
#'
#' out <- bison_solr(scientificName='Helianthus annuus', rows=800)
#' bisonmap(out)
#'
#' # Using additional solr fields
#' ## Faceting
#' bison_solr(scientificName='Helianthus annuus', rows=0, facet='true',
#'  facet.field='computedStateFips')
#'
#' ## Highlighting
#' bison_solr(scientificName='Helianthus annuus', rows=10, hl='true',
#'  hl.fl='scientificName')
#'
#' ## Use of hierarchy_homonym_string
#' bison_solr(hierarchy_homonym_string = '-202423-914154-914156-158852-')
#' ## -- This is a bit unwieldy, but you can find this string in the output of a call, like this
#' (string <- bison_solr(scientificName='Ursus americanus', rows=1)$points$hierarchy_homonym_string)
#' bison_solr(hierarchy_homonym_string = string)
#'
#' # The pointPath parameter
#' bison_solr(pointPath = '/-110.0,45.0/specimen')
#'
#' # Curl options
#' library("httr")
#' bison_solr(scientificName='Ursus americanus', callopts=verbose())
#' }

bison_solr <- function(decimalLatitude=NULL, decimalLongitude=NULL, year=NULL, providerID=NULL,
  resourceID=NULL, pointPath=NULL, basisOfRecord=NULL, eventDate=NULL, computedCountyFips=NULL,
  computedStateFips=NULL, scientificName=NULL, hierarchy_homonym_string=NULL, TSNs=NULL,
  recordedBy=NULL, occurrenceID=NULL, catalogNumber=NULL, ITIScommonName=NULL,
  kingdom=NULL, callopts=list(), verbose=TRUE, ...)
{
  url <- "http://bison.usgs.ornl.gov/solrstaging/occurrences/select/"
  qu <- bs_compact(list(decimalLatitude=decimalLatitude,
                     decimalLongitude=decimalLongitude,
                     year=year,
                     pointPath=pointPath,
                     providerID=providerID,
                     resourceID=resourceID,
                     basisOfRecord=basisOfRecord,
                     eventDate=eventDate,
                     computedCountyFips=computedCountyFips,
                     computedStateFips=computedStateFips,
                     scientificName=scientificName,
                     hierarchy_homonym_string=hierarchy_homonym_string,
                     TSNs=TSNs,
                     recordedBy=recordedBy,
                     occurrenceID=occurrenceID,
                     catalogNumber=catalogNumber,
                     ITIScommonName=ITIScommonName,
                     kingdom=kingdom))

  stuff <- list()
  for(i in seq_along(qu)){
     stuff[i] <- paste0(names(qu)[[i]],':"', qu[[i]], '"')
  }
  stuff <- paste0(stuff,collapse="+")

  args <- bs_compact(list(q=stuff, wt="json", ...))

  tt <- GET(url, query=args, c(config(followlocation=1), callopts))
  mssg(verbose, tt$url)
  stop_for_status(tt)
  out <- content(tt, as="text")

  temp <- list(
    num_found = fromJSON(out)$response$numFound,
    points = solr_parse_search(input = out, parsetype = "df"),
    highlight = solr_parse_highlight(out),
    facets = solr_parse_facets(out)
  )
  class(temp) <- "bison_solr"
  return( temp )
}

solr_parse_facets <- function(input, parsetype=NULL, concat=',')
{
  input <- fromJSON(input, simplifyVector = FALSE)

  # Facet queries
  fqdat <- input$facet_counts$facet_queries
  if(length(fqdat)==0){
    fqout <- NULL
  } else
  {
    fqout <- data.frame(term=names(fqdat), value=do.call(c, fqdat), stringsAsFactors=FALSE)
  }
  row.names(fqout) <- NULL

  # facet fields
  ffout <- lapply(input$facet_counts$facet_fields, function(x){
    data.frame(do.call(rbind, lapply(seq(1, length(x), by=2), function(y){
      x[c(y, y+1)]
    })), stringsAsFactors=FALSE)
  })

  # Facet dates
  if(length(input$facet_counts$facet_dates)==0){
    datesout <- NULL
  } else
  {
    datesout <- lapply(input$facet_counts$facet_dates, function(x){
      x <- x[!names(x) %in% c('gap','start','end')]
      x <- data.frame(date=names(x), value=do.call(c, x), stringsAsFactors=FALSE)
      row.names(x) <- NULL
      x
    })
  }

  # Facet ranges
  if(length(input$facet_counts$facet_ranges)==0){
    rangesout <- NULL
  } else
  {
    rangesout <- lapply(input$facet_counts$facet_ranges, function(x){
      x <- x[!names(x) %in% c('gap','start','end')]$counts
      data.frame(do.call(rbind, lapply(seq(1, length(x), by=2), function(y){
        x[c(y, y+1)]
      })), stringsAsFactors=FALSE)
    })
  }

  # output
  return( list(facet_queries = replacelength0(fqout),
               facet_fields = replacelength0(ffout),
               facet_dates = replacelength0(datesout),
               facet_ranges = replacelength0(rangesout)) )
}

solr_parse_highlight <- function(input, parsetype='list', concat=',')
{
  input <- fromJSON(input, simplifyVector = FALSE)
  if(parsetype=='df'){
    dat <- input$highlight
    #       highout <- data.frame(term=names(dat), value=do.call(c, dat), stringsAsFactors=FALSE)
    highout <- data.frame(cbind(names=names(dat), do.call(rbind, dat)))
    row.names(highout) <- NULL
  } else
  {
    highout <- input$highlight
  }
  return( highout )
}


solr_parse_search <- function(input, parsetype='list', concat=',')
{
  input <- fromJSON(input, FALSE)
  if(parsetype=='df'){
    dat <- input$response$docs
    dat2 <- lapply(dat, function(x){
      lapply(x, function(y){
        if(length(y) > 1 || is(y, "list")){
          paste(y, collapse=concat)
        } else { y  }
      })
    })
    datout <- do.call(rbind.fill, lapply(dat2, data.frame, stringsAsFactors=FALSE))
    datout$X_version_ <- NULL
  } else
  {
    datout <- input
  }
  return( datout )
}

# small function to replace elements of length 0 with NULL
replacelength0 <- function(x) if(length(x) < 1){ NULL } else { x }
