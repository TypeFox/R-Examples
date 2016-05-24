#' Call Asset API
#'
#' Calls NASA's Earth Imagery Assets API and returns data.frame with information on time and location of images between two dates.
#' @param key Key for API authentication.
#' @param lon Longitud of coordinate position.
#' @param lat Latitud of coordinate position.
#' @param start_date Start date to search for image. In YYYY-MM-DD format.
#' @param end_date End date to search for image. In YYYY-MM-DD format. Defaults to current system date.
#' @examples
#'\dontrun{
#' key <- "123key"
#' img <- earth_asset(key, -100.31008, 25.66779, "2016-01-01")
#'}
#' @export
earth_asset <- function(key, lon, lat, start_date, end_date = Sys.Date()){
  tryCatch({
    start_date <- as.Date(start_date)
    end_date <- as.Date(end_date)
  },
  error = function(e){
    stop("Date parameters must be YYYY-MM-DD")
  })

  # Validate a few things
  if(!is.numeric(lon)){
    stop("Lon parameter must be numeric")
  }
  if(!is.numeric(lat)){
    stop("Lat parameter must be numeric")
  }

  # fix dates
  if(Sys.Date() == end_date){
    difdate <- FALSE
  }else{
    difdate <- TRUE
  }

  h <- "https://api.nasa.gov/planetary/earth/assets?"

  query <- paste0(h,
                  "lon=", lon, "&",
                  "lat=", lat, "&",
                  "begin=", start_date, "&",
                  if(difdate){paste0("end=",end_date,
                                     "&api_key=",key)
                    }else{paste0("&api_key=",key)})
  s <- jsonlite::fromJSON(query)

  if("error" %in% names(s)){
    stop(cat(paste0("NASA API Error \n",
                    "The following is the output: ", s$error )))
  }

  # coordinates (to be used elsewhere)
  type <- "Point"
  coordinates <- paste0(as.character(lon), " ", as.character(lat))

  if(s$count<1){
    df <- data.frame("date" = "1900-01-01",
                     "id" = "NO RESULTS",
                     "type" = type,
                     "coordinates" = coordinates)
  }else{
    df <- data.frame("date" = s$results$date,
                     "id" = s$results$id,
                     "type" = type,
                     "coordinates" = coordinates)
  }
  return(df)
}
