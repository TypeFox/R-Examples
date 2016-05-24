#' OData: A helper package to access data from a web service based on OData Protocol.
#'
#' The OData package provides methods to easily call an OData web service and get data from datasets
#'
#' @docType package
#' @name OData
NULL
#> NULL

#' Get the metadata for a specific web service based on OData Protocol
#' @export
#' @title Get the metadata for a specific web service based on OData Protocol
#' @name metadata
#' @param url URL to OData ($metada can be included as parameter to get more informations)
#' @return The metadata of the OData target
#' @author BPM-Conseil
#' @examples
#' metadata("http://odata.research.microsoft.com/odata.svc/$metadata")
metadata <- function(url) {
  metadataContent <- XML::xmlParse(paste(readLines(url, warn=FALSE), collapse=""))
  root = XML::xmlRoot(metadataContent)
  return (root)
}

#' Get a list of EntitySets for a web service base on OData
#' @export
#' @title Get a list of EntitySets for a web service base on OData
#' @name entitySets
#' @param url URL to OData (needs to include $format=json as parameter)
#' @return EnitySets of the OData target
#' @author BPM-Conseil
#' @examples
#' entitySets("http://odata.research.microsoft.com/odata.svc/?$format=json")
entitySets <- function(url) {
  sets <- RJSONIO::fromJSON(paste(readLines(url, warn=FALSE), collapse=""))

  df <- data.frame(name=character(),
                   stringsAsFactors=FALSE)

  # i <- 0
  for( item in sets[2]$value){
    if (which(item != 'kind') || item[["kind"]] == "EntitySet"){
      df2 <- data.frame(name = item[1])
      df <- rbind(df, df2)
    }
  }

  return (df)
}

#' Retrieve data as list from a dataset
#' @export
#' @title Retrieve data as list from a dataset
#' @name retrieveData
#' @param url URL to dataset (needs to include $format=json as parameter)
#' @return Data from the OData web service
#' @author BPM-Conseil
#' @examples
#' retrieveData("http://odata.research.microsoft.com/odata.svc/Labs?$format=json")
retrieveData <- function(url) {
  sets <- RJSONIO::fromJSON(paste(readLines(url, warn=FALSE), collapse=""))
  return (sets)
}

#' Download a CSV file from an URL
#' @export
#' @title Download a CSV file from an URL
#' @name downloadResourceCsv
#' @param url URL pointing at a CSV file
#' @param separator Separator for the CSV file
#' @return Data from a CSV file
#' @author BPM-Conseil
downloadResourceCsv <- function(url, separator) {
  format <- "csv"
  if (format == "csv") {
    destinationFile <- "downloaded_file.csv"
  }
  else if (format == "xls") {
    destinationFile <- "downloaded_file.xls"
  }

  utils::download.file(url, destfile = destinationFile )

  csv.data <- utils::read.csv(destinationFile, sep=separator)
  return(csv.data)
}
