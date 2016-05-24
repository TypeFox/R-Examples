#' Retrieve data from the Census API
#'
#' Returns Census data for the 2010 SF1 or ACS 2013-2015 1-, 3-, and 5-Yr
#' aggregations for requested variables and geographies.
#'
#' @details See \code{vignette('censusr', package = 'censusr')} for examples.
#'
#' @param variables_to_get A character vector of the desired variable names for
#'   the Census API call, defined at \url{http://api.census.gov/}
#' @param names A character vector of the same length as \code{variables_to_get}
#'   giving the user-defined names for the variables (optional). Defaults to raw
#'   API names.
#' @param geoids A character vector of FIPS codes; must be at least to the
#'   county (5-digit) level, and can accept down to blocks (15-digit).
#' @param allgeos (optional) A string identifying the type of geography for
#'   which to collect data within the the requested \code{geoids}. Must be one
#'   of \code{c('co', 'tr', 'bg', 'bl')}. For instance, if \code{allgeos =
#'   "bg"}, will return all block groups within the given \code{geoids}.
#' @param data_source A string identifying whether the SF1 (decennial census) or
#'   ACS data is desired.
#' @param year If \code{data_source = "acs"}, the final year of the summary
#'   period. Default is \code{2013}.
#' @param period If \code{data_source = "acs"}, the length of aggregation period.
#'   Default is \code{5}, or a 5-year aggregation table.
#' @param api_key The user's Census API key (as a character string). You can get
#'   a free key from [Census](http://api.census.gov/data/key_signup.html). See
#'   \code{vignette('censusr', package = 'censusr')} to setup a default key as
#'   an environment variable.
#'
#' @return a data_frame with each requested variable at each requested geography.
#'
#' @export
call_census_api <- function(variables_to_get,
                            names = NULL,
                            geoids, allgeos = NULL,
                            data_source = c("sf1", "acs"),
                            year = 2013, period = 5,
                            api_key = NULL){

  if(Sys.getenv("CENSUS_TOKEN") == "" && is.null(api_key)){
    stop("censusr requires an API key. Request one at http://api.census.gov/data/key_signup.html")
  }
  if(is.null(api_key)) {
    api_key = Sys.getenv("CENSUS_TOKEN")
  }

  if(!is.null(allgeos)){
    if(!(allgeos %in% c("co", "tr", "bg", "bl"))){
      stop("`allgeos` must be one of c('co', 'tr', 'bg', 'bl')")
    }
  }

  # call_api_once for each requested geography
  d <- do.call(
    "rbind",
    lapply(geoids, function(geoid)
      call_api_once(variables_to_get, geoid, allgeos,
                    data_source, year, period, api_key)
    )
  )

  if(is.null(names)){
    message(
      "Returning raw variable names; pass `names` vector for user-specified names"
    )
  } else if(length(variables_to_get) != length(names)){
    warning(
      "length(names) must equal length(variables_to_get); returning raw variable names"
      )
    } else {
    names(d) <- c("geoid", names)
  }


  return(d)
}

#' Call Census API for a set of variables
#'
#' This is an internal function and is not intended for users. See instead
#' \link{call_census_api}.
#'
#' @inheritParams call_census_api
#' @param geoid A character string with a FIPS code, between 2 and 15 digits long.
#'
#' @return A code{data.frame} with the requested variables at the requested
#'   geography.
#'
#' @importFrom httr content GET
#' @importFrom dplyr select tbl_df
call_api_once <- function(variables_to_get, geoid, allgeos, data_source, year,
                          period, api_key) {

  # construct primary url depending on requested dataset
  if(data_source == "sf1"){
    # Census SF1 data
    call_start <- "http://api.census.gov/data/2010/sf1?get="
  } else if(data_source == "acs"){
    # ACS summary tables
    call_start <- paste(
      "http://api.census.gov/data/", year,
      "/acs", period, "?get=", sep = ""
    )
  }

  # construct variable url string
  var_string <- paste(variables_to_get, collapse = ",")

  # construct geo url string
  geo_string <- get_geo_url(geoid, allgeos)

  # consruct api string
  api_string = paste0("&key=", api_key)

  # assemble url
  url <- paste0(call_start, var_string, geo_string, api_string)

  # gives back a list of lists; convert to dataframe
  response <- httr::content(httr::GET(url))

  df <- data.frame(t(sapply(response, c)), stringsAsFactors = F)[-1,]
  names(df) <- response[[1]]

  # concatenate geoid
  df$geoid <- paste0(df$state,
                if("county" %in% names(df)) {df$county},
                if("tract" %in% names(df)) {df$tract},
                if("block group" %in% names(df)) {df$"block group"},
                if("block" %in% names(df)) {df$block})

  # Reorder with geoid as first column
  col_indexes <- match(variables_to_get, names(df))
  df <- dplyr::select(df, geoid, col_indexes)
  df[,-1] <- as.numeric(unlist(df[,-1]))

  return(dplyr::tbl_df(df))
}


#' Construct a geography request string from a FIPS Code
#'
#' @inheritParams call_api_once
#' @return A string with the FIPS formatted for an API request.
#'
get_geo_url <- function(geoid, allgeos) {

  split_geo <- function(geoid) {
    list(
      st = substr(geoid, 1, 2),
      co = substr(geoid, 3, 5),
      tr = substr(geoid, 6, 11),
      bg = substr(geoid, 12, 12),
      bl = substr(geoid, 12, 15)
    )
  }

  newgeo <- split_geo(geoid)
  st <- newgeo$st; co <- newgeo$co; tr <- newgeo$tr;
  bg <- newgeo$bg; bl <- newgeo$bl

  if(is.null(allgeos)) {  # if `allgeos` is not specified
    if(bl != ""){
      # blocks
      paste0(
        "&for=block:", bl,
        "&in=state:", st,
        "+county:", co,
        "+tract:", tr
      )

    } else if(bg != ""){
      # block groups
      paste0(
        "&for=block+group:", bg,
        "&in=state:", st,
        "+county:", co,
        "+tract:", tr
      )

    } else if(tr != ""){
      # tracts
      paste0(
        "&for=tract:", tr,
        "&in=state:", st,
        "+county:", co
      )
    } else {
      # counties
      paste0(
        "&for=county:", co,
        "&in=state:", st
      )
    }
  } else {  # if `allgeos` is specified
    # get `for` part
    map = data.frame(
      abbr = c("co", "tr", "bg", "bl"),
      full = c("county", "tract", "block+group", "block"),
      stringsAsFactors = F
    )
    pre = paste0("&for=",
                 map[which(map$abbr == allgeos), 'full'],
                 ":*")

    # return pre + geoid
    if(bg != ""){
      # block groups
      paste0(
        pre,
        "&in=block+group:", bg,
        "+state:", st,
        "+county:", co,
        "+tract:", tr
      )

    } else if(tr != ""){
      # tracts
      paste0(
        pre,
        "&in=tract:", tr,
        "+state:", st,
        "+county:", co
      )
    } else {
      # counties
      paste0(
        pre,
        "&in=county:", co,
        "+=state:", st
      )
    }
  }
}
