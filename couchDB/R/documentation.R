#' Functions for interacting with a couchDB database,
#' 
#'
#' This package contains various functions for interacting with couchDB  - a document database.
#' For more information about the couchDB database see: http://couchdb.apache.org.
#' @section Get connected:
#' To interact with a couchDB instance you will need to create a connection object. Communication with couchDB happens over a http protocol. A minimal connection (running on the default port and no password protection) would be:
#' 
#' \code{
#' myConn  <- couch_http_connection("localhost")
#' }
#' 
#' The variable "myConn" can now be used as parameters to other functions. 
#' 
#' For convenience a default connection can also be created with \link{couch_set_default_connection} using the same parameters.
#' 
#' Once a connection object exists, you may want to make sure it is connecting correctly with the \link{couch_ping} function. If you are properly connected the reponse should be like:
#' \preformatted{
#' Response [http://localhost:5984]
#' Status: 200
#' Content-type: text/plain; charset=utf-8
#' Size: 151 B
#' \{"couchdb":"Welcome","uuid":"c1a367c91517195b57ddafe788a72b75","version":"1.4.0","vendor"
#' \{"name":"...
#' }
#' 
#' @section Databases:
#' On a couchDB install there can be several databases (i.e. namespaces). The function \link{couch_list_databases} will provide a list of available databases on the connection provided.
#' 
#' The function \link{couch_create_database} will, similarly, allow you to create a new database or namespace on the couchDB instance.
#' @section Fetch, store and delete documents:
#' Once you have a connection and a database to work with you can fetch, store and delete documents by using the corresponding functions:
#'\itemize{
#'\item{\link{couch_fetch} and \link{couch_fetch_default}}
#'\item{\link{couch_store} and \link{couch_store_default}}
#'\item{\link{couch_delete}}
#'}
#' @author Aleksander Dietrichson
#' @docType package
#' @name couchDB
require(bitops)
require(RCurl)
require(httr)
require(rjson)