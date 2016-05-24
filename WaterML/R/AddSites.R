#' AddSites
#'
#' This function adds a table of sites to HydroServer Lite.
#' The input must be a data.frame with all required ODM site fields
#' NOTE: this only works with HydroServer Lite that implements the JSON API.
#' you must specify a valid server url, user name, and password for the HydroServer.
#' The examples here use the 'sandbox' HydroServer on http://worldwater.byu.edu/app/
#' with the username: admin and password: password.
#'
#' @import RJSONIO
#' @import httr
#' @param server The URL of the web service ending with /services or with ?wsdl,
#'  for example: http://worldwater.byu.edu/app/index.php/default/services/cuahsi_1_1.asmx?wsdl
#'  alternatively you can specify the JSON API url like:
#'  http://worldwater.byu.edu/app/index.php/default/services/api/
#' @param username The valid HydroServer Lite username, for example "admin"
#' @param password The valid HydroServer Lite password, for example "password"
#' @param sites The valid table of sites. This table must have the following 4 columns:
#' SiteCode, SiteName, Latitude, Longitude. It can also have the optional columns:
#' Elevation, SiteType, State, County, Comments.
#' @return A table of the added sites, with two extra columns:
#' SiteID (the ID assigned by the server),
#' Status (the status showing if the site was added: OK or Error). If the status is Error, then
#' the Error message with reason why the site could not be added is also shown.
#' @keywords waterml
#' @export
#' @examples
#' user <- "admin"
#' pass <- "password"
#' server <- "http://worldwater.byu.edu/app/index.php/default/services/cuahsi_1_1.asmx"
#' #make random site codes
#' random_codes = sprintf("%04d",sample(1:10000, 2))
#' random_names = paste("R","Upload", random_codes)
#' random_lats = runif(2, 35.0, 49.0) #two random latitudes inside U.S
#' random_lons = runif(2, -110.0, -70.0) #random longitudes inside U.S
#' my_sites <- data.frame(SiteCode=random_codes, SiteName=random_names,
#'                        Latitude=random_lats, Longitude=random_lons)
#'
#' added_sites <- AddSites(server, username=user, password=pass, sites=my_sites)

AddSites <- function(server, username, password, sites) {

  #check if the server is a valid url
  cuahsi <- regexpr("/cuahsi", server)
  services_api <- regexpr("/services/api", server)
  url <- NULL
  if (cuahsi > 0) {
    baseurl <- substr(server, 1, cuahsi)
    url <- paste(baseurl, "api/sites",sep="")
  } else if (services_api > 0) {
    baseurl <- substr(server, 1, services_api)
    url <- paste(baseurl, "/services/api/sites")
  } else {
    stop("The server url must contain cuahsi_1_1.asmx or ?wsdl or /services/api ")
  }

  #check if sites table has all required columns
  cols <- names(sites)
  cols.required <- c("SiteCode", "SiteName", "Latitude", "Longitude")
  cols.matched <- match(cols.required, cols)
  if (length(cols.required[is.na(cols.matched)]) > 0) {
    cols.missing <- cols.required[is.na(cols.matched)]
    msg <- paste("sites table has missing columns:", cols.missing)
    stop(msg)
  }

  i <- 1
  N <- nrow(sites)
  added.sites <- data.frame(SiteID=character(),
                            SiteName=character(),
                            SiteCode=character(),
                            Latitude=numeric(),
                            Longitude=numeric(),
                            Status=character(),
                            stringsAsFactors=FALSE)

  for(i in 1:N) {
    x <- list(
      user = username,
      password = password,
      SourceID = 1,
      SiteCode = as.character(sites$SiteCode[i]),
      SiteName = as.character(sites$SiteName[i]),
      Latitude = as.numeric(sites$Latitude[i]),
      Longitude = as.numeric(sites$Longitude[i]),
      SiteType = "Atmosphere",
      Elevation_m = 0,
      State = "a",
      County = "a",
      Comments = "a"
    )
    post.body <- RJSONIO::toJSON(x)
    print(post.body)
    #post data to server:
    response <- POST(url,
                     body = post.body,
                     add_headers("Content-Type" = "application/json")
    )
    status.code <- http_status(response)$category
    if (tolower(status.code) == "success") {
      response_status = content(response, type="application/json")
      print(response_status)
      if (response_status$status == "200 OK") {
        new.id.start <- regexpr("ID=", response_status$message)
        new.id <- substr(response_status$message, new.id.start + 3, nchar(response_status$message))
        new.site <- c(new.id, as.character(x$SiteName), as.character(x$SiteCode),
                      as.numeric(x$Latitude), as.numeric(x$Longitude), status.code)
        added.sites[nrow(added.sites)+1,] <- new.site
      } else {
        print("ERROR!")
        print(response)
      }
      #http status code is other than success...
    } else {
      new.site <- c(NA, as.character(x$SiteName), as.character(x$SiteCode),
                    as.numeric(x$Latitude), as.numeric(x$Longitude), status.code)
      added.sites[nrow(added.sites)+1,] <- new.site
    }
  }
  return (added.sites)
}
