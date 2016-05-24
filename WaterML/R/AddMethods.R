#' AddMethods
#'
#' This function adds a table of methods to HydroServer Lite.
#' The input must be a data.frame with all required ODM 'method' fields
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
#' @param methods The valid table of methods. This table must have the following columns:
#' MethodDescription, MethodLink.
#' @return A table of the added methods, with two extra columns:
#' MethodID (the ID assigned by the server),
#' Status (the status showing if the method was added: OK or Error). If the status is Error, then
#' the Error message with reason why the method could not be added is also shown.
#' @keywords waterml
#' @export
#' @examples
#' user <- "admin"
#' pass <- "password"
#' server <- "http://worldwater.byu.edu/app/index.php/default/services/cuahsi_1_1.asmx"
#' #make random site codes
#' random_description <- sprintf("R Test Method %04d",sample(1:10000, 1))
#' random_link <- "http://example.com"
#' my_methods <- data.frame(
#'   MethodDescription = random_description,
#'   MethodLink = random_link
#' )
#'
#' added_methods <- AddMethods(server, username=user, password=pass,
#'                                 methods=my_methods)

AddMethods <- function(server, username, password, methods) {

  #check if the server is a valid url
  cuahsi <- regexpr("/cuahsi", server)
  services_api <- regexpr("/services/api", server)
  url <- NULL
  if (cuahsi > 0) {
    baseurl <- substr(server, 1, cuahsi)
    url <- paste(baseurl, "api/methods",sep="")
  } else if (services_api > 0) {
    baseurl <- substr(server, 1, services_api)
    url <- paste(baseurl, "/services/api/methods")
  } else {
    stop("The server url must contain cuahsi_1_1.asmx or ?wsdl or /services/api ")
  }

  #check if table has all required columns
  cols <- names(methods)
  cols.required <- c("MethodDescription",
                     "MethodLink")
  cols.matched <- match(cols.required, cols)
  if (length(cols.required[is.na(cols.matched)]) > 0) {
    cols.missing <- cols.required[is.na(cols.matched)]
    msg <- paste("methods table has missing columns:", cols.missing)
    stop(msg)
  }

  i <- 1
  N <- nrow(methods)
  added.methods <- data.frame(MethodID=character(),
                            MethodDescription=character(),
                            MethodLink=character(),
                            Status=character(),
                            stringsAsFactors=FALSE)

  for(i in 1:N) {
    x <- list(
      user = username,
      password = password,
      MethodDescription = as.character(methods$MethodDescription[i]),
      MethodLink = as.character(methods$MethodLink[i]),
      VariableID=0
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
        new.method <-
          c(new.id,
            as.character(x$MethodDescription[i]),
            as.character(x$MethodLink[i]),
            as.character(status.code)
        )
        added.methods[nrow(added.methods)+1,] <- new.method
      } else {
        print(paste("ERROR!", status.code))
      }
      #http status code is other than success...
    } else {
      if (status.code == "server error") {
        print(response)
      }
      new.method <-
        c(NA,
          as.character(x$MethodDescription[i]),
          as.character(x$MethodLink[i]),
          as.character(status.code)
        )
      added.methods[nrow(added.methods)+1,] <- new.method
    }
  }
  return (added.methods)
}
