#' Search for and collect data from the USGS Bison API.
#'
#' @importFrom httr GET content warn_for_status stop_for_status config
#' @importFrom jsonlite fromJSON
#' @importFrom plyr ldply
#' @export
#'
#' @param species (character) A species name.
#' @param type (character) Type, one of scientific_name or common_name.
#' @param tsn (numeric) Specifies the TSN to query by. If you supply a tsn it doesn't make sense to
#' supply a species name as well. Example:162003.
#' @param start (numeric) Record to start at.
#' @param count (numeric) Number of records to return.
#' @param countyFips (character) Specifies the county fips code to geographically constrain
#'    the search to one county. Character must be supplied as a number starting
#'    with zero may lose the zero. Eg: "49015".
#' @param county (character) County name. As codes are a pain in the ass, you can put in the
#'    county name here instead of specifying a countyFips entry, and bison will
#'    attempt to look up the countyFips code. (character)
#' @param state (character) Specifies the state name to geographically constrain the search.
#'    Example: Tennessee.
#' @param aoi Specifies a WKT (Well-Known Text) polygon to geographically constrain the search.
#'  Eg.: c(-111.06 38.84,
#'         -110.80 39.377,
#'         -110.20 39.17,
#'         -110.20 38.90,
#'         -110.63 38.67,
#'         -111.06 38.84),
#'  which calls up the occurrences within the specified area. Check out the Wikipedia
#'  page here \url{http://en.wikipedia.org/wiki/Well-known_text} for an in depth
#'  look at the options, terminology, etc. (character)
#' @param aoibbox Specifies a four-sided bounding box to geographically constrain
#'    the search (using format: minx,miny,maxx,maxy). The coordinates are Spherical
#'    Mercator with a datum of WGS84. Example: -111.31,38.81,-110.57,39.21 (character)
#' @param params (character) String of parameters, one of providerID, resourceID, basisOfRecord,
#'    catalogNumber, year, computedStateFips, hierarchy_homonym_string, TSNs, recordedBy,
#'    occurrenceID, collectorNumber, provider, ownerInstitutionCollectionCode, eventDate,
#'    providedScientificName, scientificName, ITISscientificName, providedCommonName,
#'    ITIScommonName, kingdom, ITIStsn, centroid, higherGeographyID, computedCountyFips,
#'    providedCounty, calculatedCounty, stateProvince, calculatedState, countryCode. See examples.
#' @param what What to return?  One of 'all', 'summary', 'points', 'counties', 'states',
#'    'raw', or 'list'. All data is returned from the BISON API, but this parameter lets
#'    you select just the parts you want, and the rest is discarded before returning the
#'    result to you.
#' @param ... Further args passed on to httr::GET. See examples.
#'
#' @seealso \code{\link{bison_solr}} \code{\link{bison_tax}}
#'
#' @examples \dontrun{
#' bison(species="Bison bison", count=50, what='summary')
#' bison(species="Bison bison", count=50, what='points')
#' bison(species="Bison bison", count=50, what='counties')
#' bison(species="Bison bison", count=50, what='states')
#' bison(species="Bison bison", count=50, what='raw')
#' bison(species="Bison bison", count=50, what='list')
#'
#' out <- bison(species="Bison bison", count=50) # by default gets 10 results
#' out$summary # see summary
#' out$counties # see county data
#' out$states # see state data
#' bisonmap(out, tomap = "points")
#' bisonmap(out, tomap = "county")
#' bisonmap(out, tomap = "state")
#'
#' # Search for a common name
#' bison(species="Tufted Titmouse", type="common_name", what='summary')
#'
#' # Constrain search to a certain county, 49015 is Emery County in Utah
#' bison(species="Helianthus annuus", countyFips = "49015")
#'
#' # Constrain search to a certain county, specifying county name instead of code
#' bison(species="Helianthus annuus", county = "Los Angeles")
#' bison(species="Helianthus annuus", county = "Los")
#'
#' # Constrain search to a certain aoi, which turns out to be Emery County, Utah as well
#' bison(species="Helianthus annuus",
#'  aoi = "POLYGON((-111.06360117772908 38.84001566645886,
#'                  -110.80542246679359 39.37707771107983,
#'                  -110.20117441992392 39.17722368276862,
#'                  -110.20666758398464 38.90844075244811,
#'                  -110.63513438085685 38.67724220095734,
#'                  -111.06360117772908 38.84001566645886))")
#'
#' # Constrain search to a certain aoibbox, which, you guessed it, is also Emery Co., Utah
#' bison(species="Helianthus annuus", aoibbox = '-111.31,38.81,-110.57,39.21')
#'
#' # Taxonomic serial number
#' bison(tsn=162003)
#' ## If you don't have tsn's, search for a taxonomic serial number
#' library('taxize')
#' poa_tsn <- get_tsn('Poa annua')
#' bison(tsn=poa_tsn)
#'
#' # Curl debugging and other httr options, some of these examples aren't that useful, but
#' # are given for demonstration purposes
#' library("httr")
#' ## get curl verbose output to see what's going on with your request
#' bison(tsn=162003, count=1, what="points", config=verbose())
#' ## set a timeout so that the call stops after time x, compare 1st to 2nd call
#' # bison(tsn=162003, count=1, what="points", config=timeout(seconds=1))
#' # bison(tsn=162003, count=1, what="points", config=timeout(seconds=0.1))
#' ## set cookies
#' bison(tsn=162003, count=1, what="points", config=set_cookies(a = 1, b = 2))
#' ## set cookies
#' bison(tsn=162003, count=1, what="points", config=user_agent("rbison"))
#'
#' # Params - the params function accepts a number of search terms
#' ## Find the provider with ID 318.
#' bison(params='providerID:("318")')
#' ## Find all resources with id of '318,1902' OR '318,9151', with values separated by spaces.
#' bison(params='resourceID:("318,1902" "318,9151")')
#' ## Criterion may be combined using the semicolon (';') character, which translates to a logical
#' ## AND operator. Note that field names and values are case sensitive.
#' bison(params='providerID:("408" "432");resourceID:("14027")')
#' ## Search by basisOfRecord, for specimen types in this case
#' bison(params='basisOfRecord:(specimen)')
#' ## Search by computedStateFips, 01 for Alabama
#' bison(params='computedStateFips:01')
#' ## Search by ITIStsn
#' bison(params='ITIStsn:162003')
#' ## Search by countryCode
#' bison(params='countryCode:US')
#' ## Search by ITIScommonName
#' bison(params='ITIScommonName:"Canada goose"')
#' }

bison <- function(species=NULL, type="scientific_name", tsn=NULL, start=NULL, count=10,
  countyFips=NULL, county=NULL, state=NULL, aoi=NULL, aoibbox=NULL, params=NULL,
  what='all', ...)
{
  stopifnot(is.numeric(count))
  stopifnot(count >= 0)

  if(is.null(species)){
    type <- NULL
  }

  countyFips <- county_handler(county)

  if(!is.null(tsn)){
    itis <- 'itis'
    tsn <- as.numeric(as.character(tsn))
    stopifnot(is.numeric(tsn))
  } else { itis <- NULL }

  # check if param names are in the accepted list
  check_params(params)

  url <- "http://bison.usgs.ornl.gov/api/search.json"
  args <- bs_compact(list(species=species,type=type,itis=itis,tsn=tsn,start=start,count=count,
                       countyFips=countyFips,state=state,aoi=aoi,aoibbox=aoibbox,params=params))
  tt <- GET(url, query=args, ...)
  warn_for_status(tt)
  if(tt$status_code > 201){
    stopifnot(tt$headers$`content-type` == "text/html;charset=utf-8")
  } else {
    stopifnot(tt$headers$`content-type` == "application/json;charset=UTF-8")
  }
  if(tt$status_code > 201){
    res <- NA
  } else {
    out <- content(tt, as="text")
    json <- fromJSON(out, FALSE)
    what <- match.arg(what, choices=c("summary", "counties", "states", "points", "all", "raw", "list"))
    res <- switch(what,
                  summary=bison_data(json, "summary"),
                  all=bison_data(json, "all"),
                  counties=bison_data(json, "counties"),
                  states=bison_data(json, "states"),
                  points=bison_data(json, "points"),
                  raw=out,
                  list=json
    )
  }
  structure(res, class="bison")
}

check_params <- function(x){
  if(!is.null(x)){
    y <- strsplit(x, ";")[[1]]
    z <- vapply(y, function(b) strsplit(b, ":")[[1]][[1]], "", USE.NAMES = FALSE)
    check <- z %in% c('providerID','resourceID','basisOfRecord','catalogNumber','year','computedStateFips',
                      'hierarchy_homonym_string','TSNs','recordedBy','occurrenceID','collectorNumber',
                      'provider','ownerInstitutionCollectionCode','eventDate','providedScientificName',
                      'scientificName','ITISscientificName','providedCommonName','ITIScommonName','kingdom',
                      'ITIStsn','centroid','higherGeographyID','computedCountyFips','providedCounty',
                      'calculatedCounty','stateProvince','calculatedState','countryCode')
    if(!all(check)) stop("You used in an incorrect param field", call. = FALSE)
  }
}

bison_data <- function(input = NULL, datatype="summary")
{
  if(datatype=='summary'){
    tt <- data.frame(c(input[1], input$occurrences$legend))
    list(summary=tt, states=NULL, counties=NULL, points=NULL)
  } else if(datatype=="counties"){
    tt <- getcounties(input)
    list(summary=NULL, states=NULL, counties=tt, points=NULL)
  } else if(datatype=="states"){
    tt <- getstates(input)
    list(summary=NULL, states=tt, counties=NULL, points=NULL)
  } else if(datatype=="points"){
    tt <- getpoints(input)
    list(summary=NULL, states=NULL, counties=NULL, points=tt)
  } else if(datatype=="all"){
    summary=data.frame(c(input[1], input$occurrences$legend))
    counties=getcounties(input)
    states=getstates(input)
    points=getpoints(input)
    list(summary=summary, states=states, counties=counties, points=points)
  }
}

getcounties <- function(x){
  tryx <- tryCatch(x$counties$total, error = function(e) e)
  if(is(tryx, "simpleError") || is.null(tryx)){
    NULL
  } else {
    if(x$counties$total == 0){
      NULL
    } else {
      if(class(x$counties$data[[1]])=="character"){
        df <- ldply(x$counties$data)
      } else
      {
        df <- ldply(x$counties$data, function(y) data.frame(y))
      }
      names(df)[c(1,3)] <- c("record_id","county_name")
      return(df)
    }
  }
}

getstates <- function(x){
  tryx <- tryCatch(x$states$total, error = function(e) e)
  if(is(tryx, "simpleError") || is.null(tryx)){
    NULL
  } else {
    if(x$states$total == 0){
      NULL
    } else {
      df <- ldply(x$states$data, function(y) data.frame(y))
      names(df)[c(1,3)] <- c("record_id","county_fips")
      return(df)
    }
  }
}

getpoints <- function(x){
  tryx <- tryCatch(x$data, error = function(e) e)
  if(is(tryx, "simpleError")){
    NULL
  } else if(length(x$data) == 0){
      NULL
  } else {
    withlatlong <- x$data[sapply(x$data, length, USE.NAMES=FALSE) == 8]
    data_out <- ldply(withlatlong, function(y){
      y[sapply(y, is.null)] <- NA
      data.frame(y[c("name","decimalLongitude","decimalLatitude","occurrenceID",
                     "provider","basis","common_name","geo")], stringsAsFactors=FALSE)
    })
    data_out$decimalLongitude <- as.numeric(as.character(data_out$decimalLongitude))
    data_out$decimalLatitude <- as.numeric(as.character(data_out$decimalLatitude))
    return(data_out)
  }
}

county_handler <- function(x){
  if(!is.null(x)){
    numbs <- fips[grep(x, fips$county),]
    if(nrow(numbs) > 1){
      message("\n\n")
      print(numbs)
      message("\nMore than one matching county found '", x, "'!\nEnter row number of county you want (other inputs will return 'NA'):\n") # prompt
      take <- scan(n = 1, quiet = TRUE, what = 'raw')

      if(length(take) == 0)
        take <- 'notake'
      if(take %in% seq_len(nrow(numbs))){
        take <- as.numeric(take)
        message("Input accepted, took county '", as.character(numbs[take, "county"]), "'.\n")
        countyFips <- paste0(numbs[take, c("fips_state","fips_county")],collapse="")
      } else {
        countyFips <- NA
        message("\nReturned 'NA'!\n\n")
      }
    } else
      if(nrow(numbs) == 1){
        countyFips <- paste0(numbs[, c("fips_state","fips_county")],collapse="")
      } else
      { stop("a problem occurred finding the countyFips...") }
  } else { countyFips <- NULL }
  return( countyFips )
}
