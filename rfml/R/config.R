#' Set up a MarkLogic database for use with rfml.
#'
#' The function installs \href{http://docs.marklogic.com/guide/rest-dev/extensions}{REST extensions} and
#' modules needed to use the package againts a MarkLogic Server database. The function needs to be executed once
#' for each database that is going to be used with rfml. It also creates a document, /rfml/rfmlInfo.json, that
#' stores the version of the rfml package and the date the database are initiated.
#'
#' The database must have a \href{http://docs.marklogic.com/guide/admin/http}{REST server}
#' and a \href{http://docs.marklogic.com/guide/admin/databases#id_38484}{module database}.
#' It also adds a document, /rfml/rfmlInfo.json, that stores the version of the rfml
#' package and the date the database are initiated.
#'
#' The user that is used for the function need to have the rest-admin role,
#' or at least the following privileges:
#' \itemize{
#'  \item http://marklogic.com/xdmp/privileges/rest-admin
#'  \item http://marklogic.com/xdmp/privileges/rest-writer
#'  \item http://marklogic.com/xdmp/privileges/rest-reader
#'  }
#'
#' @param host The hostname or ip-adress of the MarkLogic http server. Default to localhost.
#' @param port The port number of the MarkLogic http server. 8000 is used default
#' @param adminuser The username of a user that have rights to install options. admin is default.
#' @param password The password admin is default.
#' @return Nothing if success or raise a error.
#' @examples
#' \dontrun{
#' ml.init.database("localhost", "8000", "admin", "admin")
#' }
#' @export
ml.init.database <- function(host = "localhost", port = "8000", adminuser = "admin", password = "admin") {
  # general URL, used as basis for all
  mlHost <- paste("http://", host, ":", port, sep="")
  # name of libs used
  mlLibs <- .rfmlEnv$mlLibs
  # name of extensions used
  mlExts <- .rfmlEnv$mlExts

 for (i in 1:length(mlLibs)) {
  # install the needed module library
    if (suppressWarnings(.insert.lib(mlHost, adminuser, password, mlLibs[i]))) {
      message(paste("Library ",mlLibs[i] ," is now installed on ", host, ":", port, sep=""))
    }
  }

  for (i in 1:length(mlExts)) {
    if (suppressWarnings(.insert.ext(mlHost, adminuser, password, mlExts[i]))) {
      message(paste("REST extension ",mlExts[i], " is now installed on ", host, ":", port, sep=""))
    }
  }
  suppressWarnings(closeAllConnections())
  # need to store the current version of rfml in a document in the database that is used for checks...
  rfmlVer <- as.character(packageVersion("rfml"))
  initDate <- as.character(Sys.Date())
  # Need to post a config document that can used to verify against
  mlURL <- paste(mlHost, "/v1/resources/rfml.check", sep="")
  queryArgs <- list('rs:rfmlVersion'=rfmlVer, 'rs:initDate'=initDate)
  response <- PUT(mlURL,query = queryArgs, authenticate(adminuser, password, type="digest"), accept_json())
  status_code <- response$status_code
  if (status_code != 204) {
    rContent <- content(response)
    errorMsg <- paste("statusCode: ",
                      rContent$errorResponse$statusCode,
                      ", status: ", rContent$errorResponse$status,
                      ", message: ", rContent$errorResponse$message, sep="")
    stop(paste("Ops, something went wrong.", errorMsg))
  }

  message(paste(host, ":", port, " is now ready for use with rfml",sep=""))
}

#' Remove all rfml internal files in a MarkLogic database.
#'
#' The function removes the \href{http://docs.marklogic.com/guide/rest-dev/extensions}{REST extensions} and
#' modules added with the \link{ml.init.database} function. It also removes the document, /rfml/rfmlInfo.json,
#' that stores the version of the rfml package and the date the database are initiated.
#'
#' The user that is used for the login must have the  rest-admin role, or the following privileges:
#' \itemize{
#'  \item http://marklogic.com/xdmp/privileges/rest-admin
#'  \item http://marklogic.com/xdmp/privileges/rest-writer
#'  \item http://marklogic.com/xdmp/privileges/rest-reader
#'  }
#'
#' @param host The hostname or ipadress of the MarkLogic http server. Default to localhost.
#' @param port The port number of the MarkLogic http server. 8000 is used default
#' @param adminuser The username of a user that have rights to install options. admin is default.
#' @param password The password admin is default.
#' @return Nothing if success otherwise it will raise an error.
#' @examples
#' \dontrun{
#' ml.clear.database("localhost", "8000", "admin", "admin")
#' }
#' @export
ml.clear.database <- function(host = "localhost", port = "8000", adminuser = "admin", password = "admin") {
  # general URL, used as basis for all
  mlHost <- paste("http://", host, ":", port, sep="")

  # name of libs used
  mlLibs <- .rfmlEnv$mlLibs
  # name of exts used
  mlExts <- .rfmlEnv$mlExts
  # remove the document with rfml info first while we still have the exstention left...
  mlURL <- paste(mlHost, "/v1/resources/rfml.check", sep="")
  response <- DELETE(mlURL, authenticate(adminuser, password, type="digest"), accept_json())

  for (i in 1:length(mlLibs)) {
    # install the needed module library
    if (.remove.lib(mlHost, adminuser, password, mlLibs[i])) {
      message(paste("Library ",  mlLibs[i]," is removed from ", host, ":", port, sep=""))
    }
  }

  for (i in 1:length(mlExts)) {
    # install the needed module library
    if (.remove.ext(mlHost, adminuser, password, mlExts[i])) {
      message(paste("REST extension ",  mlExts[i]," is removed from ", host, ":", port, sep=""))
    }
  }

  message(paste(host, ":", port, " cleared of rfml specific files",sep=""))
}
#' Load sample data set into MarkLogic server
#'
#' The function uploads a sample data set to MarkLogic Server and returns a ml.data.frame object.
#' Provided data sets are:
#' \itemize{
#'  \item "baskets" - sample order documents that can be used with the \link{ml.arules} function.
#'  }
#' To remove the sample use the \link{rm.ml.data.frame} on the returned ml.data.frame object.
#'
#' @param conn A \link{ml.conn-class} with a connection to a MarkLoic server
#' @param dataSet Which dataset to upload, "baskets"
#' @param name The name of the object. The data will be added to a collection with that name. If not provided the dataSet name is used.
#' @return A \link{ml.data.frame} object pointing to the uploaded dataset.
#' @examples
#' \dontrun{
#'  locConn <- ml.connect()
#'  mlBaskets <- ml.load.sample.data(locConn, "baskets")
#' }
#' @export
ml.load.sample.data <- function(conn, dataSet = "baskets", name = "") {
  if (class(conn) != "ml.conn" || missing(conn)) {
    stop("Need a valid ml.conn object. Use ml.connect to create one.")
  }
  if (dataSet == "baskets") {
    dataFolder <- "extdata/baskets"
    collection <- "baskets"
  } else {
    stop("Unknown data set!")
  }
  if (nchar(name) > 0) {
    collection <- name
  }
  suppressWarnings(rfmlCollection <- .insert.ml.data(conn, dataFolder, collection, "json", ""))
  return(ml.data.frame(conn, collection=c(rfmlCollection)));
}
#' Creates or updates a Range element index.
#'
#' The function creates or updates a \href{http://docs.marklogic.com/guide/concepts/indexing#id_51573}{range element index}
#' on the underlying element/property of a \link{ml.data.frame} field.
#' The user that is used for the login needs the manage-admin role, or the following privilege:
#' \itemize{
#'  \item http://marklogic.com/xdmp/privileges/manage-admin
#'  }
#'
#' The function only creates and updates range index on a XML element or JSON property based on the \link{ml.data.frame} field.
#' Information about the field can be shown by \code{mlDataFrame$itemField}, where mlDataFrame is a \link{ml.data.frame} object
#' and itemField is the name of the field. Indexes created with this function will always have range-value-positions equal true.
#'
#' @param x a ml.data.frame field that the index will be created on
#' @param scalarType An atomic type specification. "string" is default
#' @param collation For scalarType = string, you can use a different collation than the default. Default is "http://marklogic.com/collation/"
#' @param namespaceUri The namespace URI of the XML element, if JSON ignore. Default is empty.
#' @param database The name of the database to create the index in. "Documents" is default.
#' @param host The hostname or ipadress of the MarkLogic Manage server. Default is the same as used for conn
#' @param port The port number of the MarkLogic Manage server. 8002 is used default
#' @param adminuser The username of a user that have rights to create index. Default is the same as used for conn
#' @param password The password. Default is the same as used for conn.
#' @param conn A \link{ml.conn-class} with a connection to a MarkLoic server. Optional.
#' @return The function will raise a error if something goes wrong.
#' @export
ml.add.index <- function(x, scalarType= "string", collation = "http://marklogic.com/collation/",
                         namespaceUri = "",database = "Documents", host = "", port = "8002",
                         adminuser = "", password = "", conn = NA) {

  if (class(conn) != "ml.conn" || missing(conn)) {
    hasConn <- FALSE
  } else {
    hasConn <- TRUE
  }

  # Need to check that either database etc is set or there is a valid conn object...
  if (!is.ml.col.def(x)) {
    stop("x needs to be a ml.col.def object!")
  }
  # get connection imformation
  if (nchar(adminuser) > 0) {
    username <- adminuser
  } else {
    if (!hasConn) {
      stop("You need to provide a username value or a valid ml.conn object for the conn parameter")
    } else {
      username <- conn@.username
    }
  }
  if (nchar(password) > 0) {
    pwd <- password
  } else {
    if (!hasConn) {
      stop("You need to provide a password value or a valid ml.conn object for the conn parameter")
    } else {
      key <- .rfmlEnv$key[[conn@.id]]
      pwd <- rawToChar(PKI::PKI.decrypt(conn@.password, key))
    }
  }
  if (nchar(host) > 0) {
    mlhost <- host
  } else {
    if (!hasConn) {
      stop("You need to provide a host value or a valid ml.conn object for the conn parameter")
    } else {
      mlhost <- conn@.host
    }

  }

  # general URL, used as basis for all
  mlHost <- paste("http://", mlhost, ":", port, sep="")
  mlURL <- paste(mlHost, "/manage/v2/databases/", database, "/properties", sep="")

  # first we need to get the existing properties because we need existing range-element-index, otherwise
  # we will replace them
  response <- GET(mlURL, authenticate(username, pwd, type="digest"), accept_json())
  rContent <- content(response) #, as = "text"
  #
  localname <- x@.org_name
  indexJson <- ''
  if (length(rContent$`range-element-index`) > 0 ) {
    # build up the properties for existing range-element indexes
    exRange <- rContent$`range-element-index`
    indexJsonStart <- '{"range-element-index": ['
    for (i in 1:length(exRange)) {
      # if the index we are creating already exist we will replace it
      if (exRange[[i]]$localname == localname) {
        next()
      }
      if (nchar(indexJson) > 0) {
        indexJson <- paste(indexJson, ",", sep="")
      }
      # {"range-element-index": [{"scalar-type": "string", "namespace-uri":"","localname":"productName","collation":"http://marklogic.com/collation/","range-value-positions":true,"invalid-values":"reject"},{"scalar-type": "dateTime", "namespace-uri":"","collation":"","localname":"date","range-value-positions":true,"invalid-values":"reject"}]}
      indexJson <- paste(indexJson, '{"scalar-type": "', exRange[[i]]$`scalar-type`,
                         '", "namespace-uri":"',exRange[[i]]$`namespace-uri` ,'","localname":"',exRange[[i]]$localname,
                         '","collation": "', exRange[[i]]$collation ,
                         '","range-value-positions":', tolower(exRange[[i]]$`range-value-positions`),
                         ', "invalid-values":"',exRange[[i]]$`invalid-values`,'"}', sep="")


    }
    indexJsonEnd <- ']}'
  } else {
    indexJsonStart <- '{"range-element-indexes": {"range-element-index":'
    indexJsonEnd <- '}}'
  }
  if (nchar(indexJson) > 0) {
    indexJson <- paste(indexJson, ",", sep="")
  }
   indexJson <- paste(indexJson, '{"scalar-type": "', scalarType,
                      '", "namespace-uri":"',namespaceUri ,'","localname":"', localname,'"' ,  sep="")

   # Only string that uses collation, but the collation attribute needs to be provided
   if (scalarType != "string") {
     collation <- ''
   }
   indexJson <- paste( indexJson, ',"collation": "', collation , '"' ,  sep="")
  indexJson <- paste( indexJson,',"range-value-positions":true}', sep="")
  indexJson <- paste(indexJsonStart, indexJson, indexJsonEnd, sep="")
  response <- PUT(mlURL, authenticate(username, pwd, type="digest"), body=indexJson, encode = "json", content_type_json(),accept_json())
  if(response$status_code != 204) {
    rContent <- content(response, as = "text")
    errorMsg <- paste("return message: ",
                      rContent, sep="")
    stop(paste("Ops, something went wrong.", errorMsg, "\n body: ", indexJson))
  }
  message(paste("Range element index created on ", localname,sep=""))
}
