
#' Attempt to geocode a string
#'
#' Enter a string and if found, the latitude and longitude is returned using the HERE API
#' @param search A string to search
#' @param App_id App_id to use the production HERE API. Get one here... http://developer.here.com/get-started. If left blank, will default to demo key with an unknown usage limit.
#' @param App_code App_code to use the production HERE API. Get one here... http://developer.here.com/get-started. If left blank, will default to demo key with an unknown usage limit.
#' @return A list containing Latitude and Longitude if found, NA otherwise
#' @keywords geocode
#' @export
#' @examples
#' geocodeHERE_simple("chicago")
#' geocodeHERE_simple("wrigley field chicago IL")
#' geocodeHERE_simple("233 S Wacker Dr, Chicago, IL 60606")
geocodeHERE_simple <- function(search, App_id="", App_code=""){
  if(!is.character(search)){stop("'search' must be a character string")}
  if(!is.character(App_id)){stop("'App_id' must be a character string")}
  if(!is.character(App_code)){stop("'App_code' must be a character string")}

  if(App_id=="" & App_code==""){    
    App_id <- "inwresuveWra5ebewaSweh"
    App_code <- "zBWCuMTr-PrXwr6pc5uqLg"
    base_url <- "http://geocoder.cit.api.here.com/6.2/geocode.json"
  }else{
    base_url <- "http://geocoder.api.here.com/6.2/geocode.json"
  }

  a <- httr::GET(base_url, query=list(app_id = App_id,
                                      app_code = App_code,
                                      searchtext = search))
  response <- httr::content(a)

  if(length(response$Response$View) > 0){
    ret <- response$Response$View[[1]]$Result[[1]]$Location$DisplayPosition
  }else{
    ret <- NA
  }
  return(ret)
}


