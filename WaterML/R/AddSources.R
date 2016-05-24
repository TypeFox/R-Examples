#' AddSources
#'
#' This function adds a table of sources to HydroServer Lite.
#' The input must be a data.frame with all required ODM 'Source' fields
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
#' @param sources The valid table of sources. This table must have the following columns:
#' Organization, Description, SourceLink, ContactName, ContactPhone, ContactEmail,
#' Address, City, State, Zipcode, Citation, MetadataID.
#' @return A table of the added sources, with two extra columns:
#' SourceID (the ID assigned by the server),
#' Status (the status showing if the source was added: OK or Error). If the status is Error, then
#' the Error message with reason why the source could not be added is also shown.
#' @keywords waterml
#' @export
#' @examples
#' user <- "admin"
#' pass <- "password"
#' server <- "http://worldwater.byu.edu/app/index.php/default/services/cuahsi_1_1.asmx"
#' #make random source codes
#' random_id = sample(1:10000,size=1)
#' random_name = paste("R test source", random_id)
#' my_sources <- data.frame(
#'   Organization = random_name,
#'   Description = paste("Uploaded from R:",random_name),
#'   SourceLink = paste("http://", random_id, sep=""),
#'   ContactName = random_name,
#'   ContactPhone = "012-345-6789",
#'   ContactEmail = "test<at>gmail.com",
#'   Address = random_name,
#'   City = random_name,
#'   State = random_name,
#'   Zipcode = random_id * 10,
#'   Citation = paste("Uploaded from R as a test:", random_name),
#'   MetadataID = 10
#' )
#'
#' added_sources <- AddSources(server, username=user, password=pass,
#'                                 sources=my_sources)

AddSources <- function(server, username, password, sources) {

  #check if the server is a valid url
  cuahsi <- regexpr("/cuahsi", server)
  services_api <- regexpr("/services/api", server)
  url <- NULL
  if (cuahsi > 0) {
    baseurl <- substr(server, 1, cuahsi)
    url <- paste(baseurl, "api/sources",sep="")
  } else if (services_api > 0) {
    baseurl <- substr(server, 1, services_api)
    url <- paste(baseurl, "/services/api/sources")
  } else {
    stop("The server url must contain cuahsi_1_1.asmx or ?wsdl or /services/api ")
  }

  #check if table has all required columns
  cols <- names(sources)
  cols.required <- c("Organization",
                     "Description",
                     "SourceLink",
                     "ContactName",
                     "ContactPhone",
                     "ContactEmail",
                     "Address",
                     "City",
                     "State",
                     "Zipcode",
                     "Citation",
                     "MetadataID"
                    )
  cols.matched <- match(cols.required, cols)
  if (length(cols.required[is.na(cols.matched)]) > 0) {
    cols.missing <- cols.required[is.na(cols.matched)]
    msg <- paste("sources table has missing columns:", cols.missing)
    stop(msg)
  }

  i <- 1
  N <- nrow(sources)
  added.sources <- data.frame(
    SourceID=character(),
    Organization=character(),
    Description=character(),
    SourceLink=character(),
    ContactName=character(),
    ContactPhone=character(),
    ContactEmail=character(),
    Address=character(),
    City=character(),
    State=character(),
    Zipcode=character(),
    Citation=character(),
    MetadataID=numeric(),
    Status=character(),
    stringsAsFactors=FALSE
  )

  for(i in 1:N) {
    x <- list(
      user = username,
      password = password,
      organization=as.character(sources$Organization[i]),
      description=as.character(sources$Description[i]),
      link=as.character(sources$SourceLink[i]),
      name=as.character(sources$ContactName[i]),
      phone=as.character(sources$ContactPhone[i]),
      email=as.character(sources$ContactEmail[i]),
      address=as.character(sources$Address[i]),
      city=as.character(sources$City[i]),
      state=as.character(sources$State[i]),
      zipcode=as.character(sources$Zipcode[i]),
      citation=as.character(sources$Citation[i]),
      metadata=as.numeric(sources$MetadataID[i])
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
        new.source <-
          c(new.id,
            as.character(x$organization[i]),
            as.character(x$description[i]),
            as.character(x$link[i]),
            as.character(x$name[i]),
            as.character(x$phone[i]),
            as.character(x$email[i]),
            as.character(x$address[i]),
            as.character(x$city[i]),
            as.character(x$state[i]),
            as.character(x$zipcode[i]),
            as.character(x$citation[i]),
            as.numeric(x$metadata[i]),
            status.code
        )
        added.sources[nrow(added.sources)+1,] <- new.source
      } else {
        print(paste("ERROR!", status.code))
      }
      #http status code is other than success...
    } else {
      if (status.code == "server error") {
        print(response)
      }
      new.source <-
        c(NA,
          as.character(x$organization[i]),
          as.character(x$description[i]),
          as.character(x$link[i]),
          as.character(x$name[i]),
          as.character(x$phone[i]),
          as.character(x$email[i]),
          as.character(x$address[i]),
          as.character(x$city[i]),
          as.character(x$state[i]),
          as.character(x$zipcode[i]),
          as.character(x$citation[i]),
          as.numeric(x$metadata[i]),
          status.code
        )
      added.sources[nrow(added.sources)+1,] <- new.source
    }
  }
  return (added.sources)
}
