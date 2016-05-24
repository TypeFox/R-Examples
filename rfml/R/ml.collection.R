#' Lists all collections in a MarkLogic Database.
#'
#'
#' @param conn A \link{ml.conn-class} object created by \link{ml.connect}
#' @param query Limit the collections based on a query. For more information about syntax see \link{ml.data.frame}
#' @examples
#' \dontrun{
#'  library(rfml)
#'  localConn <- ml.connect()
#'  ml.collections(localConn)
#'  }
#' @export
ml.collections <- function (conn, query="")
{
  if (class(conn) != "ml.conn" || missing(conn)) {
    stop("Need a valid ml.conn object. Use ml.connect to create one.")
  }
  # get data from ML
  # need to check that the key exits...
  key <- .rfmlEnv$key[[conn@.id]]
  password <- rawToChar(PKI::PKI.decrypt(conn@.password, key))
  username <- conn@.username

  mlHost <- paste("http://", conn@.host, ":", conn@.port, sep="")
  mlSearchURL <- paste(mlHost, "/v1/resources/rfml.collection", sep="")


  # These are the arguments that are common to all calls to MarkLogic
  queryArgs <- list('rs:q'=query)

  # do a search
  response <- GET(mlSearchURL, query = queryArgs, authenticate(username, password, type="digest"), accept_json())

  # get the content
  rContent <- content(response)

  if(response$status_code != 200) {
    errorMsg <- paste("statusCode: ",
                      rContent$errorResponse$statusCode,
                      ", status: ", rContent$errorResponse$status,
                      ", message: ", rContent$errorResponse$message, sep="")
    stop(paste("Ops, something went wrong.", errorMsg))
  }
  return(rContent)


}

#' Retrives information about a collection
#'
#' The function extracts the structure of the documents belonging to a collection
#' based on a sample it also estimates the number of documents that belongs to the collection.
#'
#' @param conn A \link{ml.conn-class} object created by \link{ml.connect}
#' @param collection A string woth the name of the collection
#' @examples
#' \dontrun{
#'  library(rfml)
#'  localConn <- ml.connect()
#'  ml.collection.info(localConn, "iris")
#'  }
#' @export
ml.collection.info <- function (conn,collection)
{
  if (class(conn) != "ml.conn" || missing(conn)) {
    stop("Need a valid ml.conn object. Use ml.connect to create one.")
  }
  if (missing(collection)) {
    stop("Need to provide a collection name")
  }
  # get data from ML
  # need to check that the key exits...
  key <- .rfmlEnv$key[[conn@.id]]
  password <- rawToChar(PKI::PKI.decrypt(conn@.password, key))
  username <- conn@.username

  mlHost <- paste("http://", conn@.host, ":", conn@.port, sep="")
  mlSearchURL <- paste(mlHost, "/v1/resources/rfml.collection", sep="")


  # These are the arguments that are common to all calls to MarkLogic
  queryArgs <- list('rs:collection'=collection)

  # do a search
  response <- GET(mlSearchURL, query = queryArgs, authenticate(username, password, type="digest"), accept_json())

  # get the content
  rContent <- content(response)

  if(response$status_code != 200) {
    errorMsg <- paste("statusCode: ",
                      rContent$errorResponse$statusCode,
                      ", status: ", rContent$errorResponse$status,
                      ", message: ", rContent$errorResponse$message, sep="")
    stop(paste("Ops, something went wrong.", errorMsg))
  }
  # return(rContent)
  nrows <- as.integer(rContent$nrows)
  fieldList <- rContent$dataFrameFields
  fieldNames <- c()
  fieldTypes <- c()
  fieldOrgNames <- c()
  fieldOrgXPaths <- c()
  fieldFormat <- c()
  fieldXmlns <- c()
  for (i in 1:length(fieldList)) {
    fieldNames[i] <-  as.character(attributes(fieldList[i]))
    fieldTypes[i] <- fieldList[[i]]$fieldType
    fieldOrgNames[i] <- fieldList[[i]]$orgField
    fieldOrgXPaths[i] <- fieldList[[i]]$orgPath
    fieldFormat[i] <- fieldList[[i]]$orgFormat
    if (!is.null(fieldList[[i]]$xmlns)) {
      fieldXmlns[i] <- fieldList[[i]]$xmlns
    } else {
      fieldXmlns[i] <- ""
    }
  }
  fields <- data.frame(name = fieldOrgNames, xpath=fieldOrgXPaths, docFormat = fieldFormat, xmlns = fieldXmlns)
  collInfo <- list(name=collection, nrows=nrows, fields=fields)
  return(collInfo)
  # Name
}
