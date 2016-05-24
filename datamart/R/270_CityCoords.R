#' Longitude and Latitude for Cities
#'
#' The UrlData object provides a resource 'CityCoordinates' that
#' takes a \code{city} parameter and returns a two-valued vector with
#' latitude and longitude.
#'
#' Based on Stackoverflow solution by Jochem Donkers.
#'
#' @param country   two-character country code, default 'DE'
#'
#' @return an object of class UrlData
#' @references 
#' \href{http://stackoverflow.com/questions/13905098/how-to-get-the-longitude-and-latitude-coordinates-from-a-city-name-and-country-i}{Stackoverflow}
#' @export
city_coords <- function(country="DE") urldata(
  resource="CityCoordinates",
  template="http://nominatim.openstreetmap.org/search?city=$(city)&countrycodes=$(country)&limit=9&format=json",
  country=country,
  city=NA,
  extract.fct=RJSONIO::fromJSON,
  transform.fct=function(x) if(is.vector(x)) c(lat=as.numeric(x[[1]]$lat), lon=as.numeric(x[[1]]$lon)) else c(NA, NA)
)


#' Reverse Geocoding: Adresse from Longitude and Latitude
#'
#' The UrlData object provides a resource 'AddressLookup' which takes two parameters
#' \code{lat} and \code{lon} and returns an approximate adress for this coordinates.
#'
#' Based on Stackoverflow solution by Jochem Donkers.
#'
#'
#' @return an object of class UrlData
#' @references 
#' \href{http://stackoverflow.com/questions/13905098/how-to-get-the-longitude-and-latitude-coordinates-from-a-city-name-and-country-i}{Stackoverflow}
#' @export
address_lookup <- function() urldata(
    resource="AddressLookup",
    template="http://nominatim.openstreetmap.org/reverse?format=json&lat=$(lat)&lon=$(lon)&zoom=18&addressdetails=1",
    lat=NA,
    lon=NA,
    extract.fct=RJSONIO::fromJSON,
    transform.fct=function(x) list(
        road=tryCatch(x$address[["road"]], error=function(e) NA), 
        city=tryCatch(x$address[["city"]], error=function(e) NA), 
        postcode=tryCatch(x$address[["postcode"]], error=function(e) NA), 
        country=tryCatch(x$address[["country"]], error=function(e) NA), 
        house_number=tryCatch(x$address[["house_number"]], error=function(e) NA)
    )
)

# address_lookup_ex <- function(verbose=TRUE) {
    # data(kw_europa_adresse)
    # osm <- address_lookup()
    # res <- transform(kw_europa_adresse, road=NA, city=NA, postcode=NA, country=NA, house_number=NA)
    # for (i in 1:nrow(kw_europa_adresse)) {
        # if(verbose) cat("looking up item ", i, " of ", nrow(kw_europa_adresse), "...\n")
        # addr_details <- query(osm, "AddressLookup", lat=kw_europa_adresse[i, "lat"], lon=kw_europa_adresse[i, "lng"])
        # res[i, c("road", "city", "postcode", "country", "house_number")] <- addr_details[c("road", "city", "postcode", "country", "house_number")]
    # }
    # return(res)
    
# }