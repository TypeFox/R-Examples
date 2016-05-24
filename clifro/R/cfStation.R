#' @include dataFrame.R
# Internals ---------------------------------------------------------------
#
# This is a function to process multiple HTTP requests concurrently to speed up
# elapsed time.
#' @importFrom RCurl getCurlHandle basicTextGatherer curlOptions curlSetOpt push 
#' complete getCurlMultiHandle
getURIs = function(uris, multiHandle = getCurlMultiHandle()) {
    content = vector("list", length(uris))
    curls = vector("list", length(uris))
    
    for(i in seq_along(uris)) {
        curl = getCurlHandle()
        content[[i]] = basicTextGatherer()
        opts = curlOptions(URL = uris[i], writefunction = content[[i]]$update)
        curlSetOpt(.opts = opts, curl = curl)
        multiHandle = push(multiHandle, curl)
    }
    
    complete(multiHandle)
    lapply(content, function(x) x$value())
}

# Create tidy names
#
# This function is used for the show method to create nicer looking row names 
# that don't take up the whole width of the page.
#
# names   : the names
# max_len : max number of characters in the name excluding spaces
tidy_names = function(names, max_len = 20){
  names = gsub("\\.+", " ", names)
  names_list = strsplit(names, " ")
  cumsum_char = lapply(names_list, function(x) cumsum(nchar(x)))
  which_lt = lapply(cumsum_char, "<", max_len)
  nice_names = sapply(mapply("[", names_list, which_lt, SIMPLIFY = FALSE),
                      paste, collapse = " ")
  cut_off = sapply(lapply(which_lt, "!"), any)
  nice_names[cut_off] = paste(nice_names[cut_off], "[..]")
  nice_names
}

# cfStation class ---------------------------------------------------------

# The Clifro Station Class
# 
# This class represents cliflo stations that can be passed into the final 
# query.
# 
#' @importFrom methods setClass
setClass("cfStation", contains = "dataFrame")

# Initializer -------------------------------------------------------------

#  Initializer for the cfStation datatype
#
# The main aim for the initializer is to order the open stations by opening date
# and the closed stations by closing date and remove duplicated stations.
# 
#' @importFrom methods setMethod
#' @importFrom lubridate now
setMethod("initialize", "cfStation", function(.Object, df){
    
    if (any(duplicated(df[3])))
        df = df[ - which(duplicated(df[3])), ]
    df_list = as.list(df)
    names(df_list) = NULL
    
    .Object@names = c("name", "network", "agent",
                      "start", "end", "open", 
                      "distance", "lat", "lon")
    
    start_diff = now() - df_list[[4]]
    end_diff = now() - df_list[[5]]
    ordered_stations = order(end_diff, -start_diff)
    .Object@.Data = lapply(df_list, "[", ordered_stations)
    .Object@row.names = paste(seq_along(df_list[[1]]))
    return(.Object)
})

# Constructor -------------------------------------------------------------

#' The Clifro Station Object
#' 
#' Create a \code{cfStation} object containing station information for one or
#' more CliFlo stations.
#' 
#' A \code{cfStation} object is created by the constructor function 
#' \code{cf_station}. The unique agent numbers of the stations are all that is 
#' required to create a \code{cfStation} object using the \code{cf_station} 
#' function. The rest of the station information including the name, network and 
#' agent ID, start and end dates, coordinates, as well as other data is scraped 
#' from CliFlo.
#' 
#' This function is used for when the agent numbers are already known. For help 
#' creating \code{cfStation} objects when the agent numbers are unknown see the
#' \code{\link{cf_find_station}} function.
#' 
#' @param ... comma separated agent numbers
#' 
#' @importFrom selectr querySelector querySelectorAll
#' @importFrom XML xmlValue htmlParse
#' @importFrom methods new
#' @importFrom lubridate dmy with_tz round_date now
#' @rdname cfStation-class
#' @name cfStation-class
#' @aliases cfStation
#' @aliases cf_station
#' @return \code{cfStation} object
#' @export
#' @seealso \code{\link{cf_find_station}} for creating \code{cfStation} objects
#' when the agent numbers are not known and \code{vignette("cfStation")} 
#' for working with clifro stations including spatial plotting in \R. For saving
#' \code{cfStation} objects as KML files refer to the vignette or 
#' \code{\link{cf_save_kml}}.
#' @examples
#' \dontrun{
#' # Create a cfStation object for the Leigh 1 and 2 Ews stations
#' leigh.st = cf_station(1339, 1340)
#' leigh.st
#' 
#' # Note, this can also be achieved using the '+' operator
#' leigh.st = cf_station(1339) + cf_station(1340)
#' leigh.st
#' 
#' # Add another column showing how long the stations have been open for
#' leigh.df = as(leigh.st, "data.frame")
#' leigh.df$ndays = with(leigh.df, round(end - start))
#' leigh.df
#' 
#' # Save the stations to the current working directory as a KML to visualise 
#' # the station locations
#' cf_save_kml(leigh.st)
#' }
cf_station = function(...){
  agent = c(...)
  if (length(agent) == 0)
    agent = 3925
  
  agent = suppressWarnings(as.numeric(agent))
  if (any(is.na(agent)))
    stop("agents must be numeric")
  
  if (any(agent %% 1 != 0))
    stop("agent numbers must be integer")
  
  uris = paste0("http://cliflo.niwa.co.nz/pls/niwp/wstn.stn_details?cAgent=", 
                agent)
  
  if (any(duplicated(agent))){
    uris = uris[!duplicated(agent)]
    agent = agent[!duplicated(agent)]
  }
  
  stations_xml = getURIs(uris)
  station_details = lapply(lapply(stations_xml, htmlParse), function(z){
    station_details_xml = querySelectorAll(z, "td.extrdata:nth-child(2)")
    if (length(station_details_xml) == 0)
      NA
    else
      sapply(station_details_xml, xmlValue)[c(1:5, 10:12)]
  })
  
  which.na = as.numeric(which(is.na(station_details)))
  if (length(which.na) == length(agent))
    stop("the agent numbers do not represent any CliFlo stations")
  
  if (length(which.na) != 0){
    if (length(which.na) == 1){
      message(paste("agent number", agent[which.na], "does not exist - ignoring"))
    } else {
      message(paste("agent numbers", paste(agent[which.na], collapse = ", "), 
                    "do not exist - ignoring"))
    }
    station_details = station_details[-which.na]
  }
  
  start_date = dmy(as.character(sapply(station_details, "[", 6)), 
                   tz = "Pacific/Auckland")
  end_date = as.character(sapply(station_details, "[", 7))
  
  open_station = end_date == "-"
  final_date = rep(with_tz(round_date(now(), "day"), 
                           tzone = "Pacific/Auckland"), length(station_details))
  
  if (any(!open_station))
    final_date[!open_station] = dmy(end_date[!open_station], 
                                    tz = "Pacific/Auckland")
  
  ## options(stringsAsFactors = FALSE)
  new("cfStation", 
      data.frame(
        name = as.character(sapply(station_details, "[", 3)),
        network = as.character(sapply(station_details, "[", 2)),
        agent = as.numeric(sapply(station_details, "[", 1)),
        start.date = start_date,
        end.date = final_date,
        open.station = open_station,
        distance = numeric(length(open_station)), 
        lat = as.numeric(sapply(station_details, "[", 4)),
        lon = as.numeric(sapply(station_details, "[", 5)),
        check.names = TRUE, stringsAsFactors = FALSE))
}

# Methods -----------------------------------------------------------------

#' @importFrom methods setMethod
setMethod("show", 
          signature(object = "cfStation"),
          function(object){
            cfstation_df = data.frame(object)
            rownames(cfstation_df) = paste0(object@row.names, ")")
            print(cfstation_df)
          })

#' @importFrom methods setMethod new
#' @rdname clifroAdd
#' @aliases +,cfStation,cfStation-method
setMethod("+",
          signature(e1 = "cfStation",
                    e2 = "cfStation"), 
          function(e1, e2){
              new.obj = new("cfStation",
                            data.frame(
                            name = c(as.character(e1@.Data[[1]]), 
                                     as.character(e2@.Data[[1]])),
                            network = c(as.character(e1@.Data[[2]]), 
                                        as.character(e2@.Data[[2]])),
                            agent = c(e1@.Data[[3]], e2@.Data[[3]]),
                            start_date = c(e1@.Data[[4]], e2@.Data[[4]]),
                            end_date = c(e1@.Data[[5]], e2@.Data[[5]]),
                            open_station = c(e1@.Data[[6]], e2@.Data[[6]]),
                            distance = c(e1@.Data[[7]], e2@.Data[[7]]),
                            latitude = c(e1@.Data[[8]], e2@.Data[[8]]),
                            longitude = c(e1@.Data[[9]], e2@.Data[[9]]),
                            stringsAsFactors = FALSE)
              )
              return(new.obj)
          })

#' @importFrom methods setMethod
#' @rdname Extract
#' @aliases [,cfStation,ANY,ANY,ANY-method
setMethod("[",
          signature(x = "cfStation"),
          function (x, i, j, drop = TRUE)
          {
            if (!missing(j))
              x = data.frame(x)[i, j, drop = drop]
            else{
              x@.Data = lapply(x@.Data, "[", i)
              x@row.names = paste(seq_along(i))
            }
            x
          }
)
