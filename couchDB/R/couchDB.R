####################################################################
# Couch R Client
# Programmer: Aleksander Dietrichson
# Purpose: interface with couchDB
# Date: 2013-11-25
require(bitops)
require(RCurl)
require(httr)
require(rjson)
#'@import bitops
#'@import RCurl
#'@import httr
#'@import rjson
#description formats name-value pair as url parameter
couch_default_database=NULL;
couch_default_connection=NULL;
url_param <- function(name, value) {
  if(is.null(value)) {
    NULL
  } else {
    if(is.logical(value)) {
      paste(name,"=",curlEscape(tolower(value)), sep="")
    } else {
      paste(name,"=",curlEscape(value), sep="")
    }
  }
}

make_query_string <- function(params) {
  params2 <- paste(params, collapse="&")
  paste("?", params2, sep="")
}

couch_base_url <- function(conn) {
  #TODO: update here for password authentication
  if(conn$secure == TRUE) {
    proto <- "https"
  } else {
    proto <- "http"
  }
  if(!is.null(conn$user)){
    auth <- paste0(conn$user,":",conn$password,"@"); 
  } else {
    auth <-""; 
  }
  base_url <- paste0(proto, "://",auth, conn$couch_http_host, ":", conn$couch_http_port)
  base_url
}

couch_fetch_url <- function(conn,database,key=NULL,opts=NULL){
  url <- paste(couch_base_url(conn), database, key, sep="/");
  url;
}

#'@title Get attachment url
#'@description Get the url for a specific attachment, This is sometimes useful for direct reads to functions, in lieu of storing tempfiles.
#'@param conn A connection object
#'@param database The database name
#'@param key Document key
#'@param attachment Name of the attachment
#'@export
couch_fetch_attachment_url <- function(conn,database,key=NULL,attachment=NULL){
  url <- paste(couch_base_url(conn), database, key,attachment, sep="/");
  url;
}

couch_store_url <- function(conn, database, key=NULL, opts=NULL) {
  if(is.null(key)) {
    url <- paste(couch_base_url(conn), database, sep="/")
  } else {
    url <- paste(couch_base_url(conn), database, key, sep="/")
    #check for revision number
    revision_number <- HEAD(url)$header$etag;
    if(!is.null(revision_number)) url <- paste0(url,"?rev=",gsub("\"","",revision_number));
  }
  url
}

couch_get_headers <- function(conn,database,key){
  path <- couch_store_url(conn,database,key);
  result <- HEAD(path)$headers;
  result;
}
couch_delete_url <- function(conn, database, key, myOpts=NULL) {
  url <- paste("databases", database, "keys", key)
  paste(url, sep="")
}
couch_stats_url <- function(conn) {
  paste(couch_base_url(conn), "/stats", sep="")
}
#TODO: update this
couch_mapred_url <- function(conn) {
  paste(couch_base_url(conn), "/mapred", sep="")
}
#'@title Set default connection
#'@description Sets up a couchDB connection to use as default
#' @param host The IP address of the couchDB instance
#' @param port The port to connect to
#' @param https Should a ssl protocol be used?
#' @export
couch_set_default_connection <- function(host,port=5984, https=FALSE){
  couch_default_connection <<- couch_http_connection(host=host,port=port,https=https);
}

#'@title Set a database as default document store
#'@description Specifies a database to write to on a couch connection by default.
#'@param database the database to use as default (String);
#'@export
couch_set_default_database <- function(database){
  if(is.character(database)){
    couch_default_database <<- database;
  } else {
    warning("Default database could not be set.");
  }
}

#title couch_list_databases_url
#description Format the url for fetching database-list
#for internal use
#param conn: A couchDB connection object.
couch_list_databases_url <- function(conn) {
  paste(couch_base_url(conn), "_all_dbs", sep="/")
}

### operations

#' @title Connection to couchDB
#' @description Creates a connection object on the host and ports provided
#' @param host The IP address of the couchDB instance
#' @param port The port to connect to
#' @param https Should a ssl protocol be used
#' @param user Username on the database server
#' @param password Password for the database server
#' @examples \dontrun{ 
#'    myConn <- couch_http_connection(host="localhost");
#' } 
#' @export
couch_http_connection <- function(host, port=5984, https=FALSE, user=NULL,password=NULL) {
  
  conn <- list(couch_http_host = host, couch_http_port = port, secure=https, user=user,password=password);
  class(conn) <- "couch_connection"
  conn
  
}


#' @title Print method for couchDB connection object
#' @description Prints the couchDB connection object.
#' @param conn a couchDB connection object
#' @method print couch_connection
print.couch_connection <- function(conn) {
  if(conn$secure== TRUE) {
    proto <- "https"
  } else {
    proto <- "http"
  }
  p <- paste(proto, conn$couch_http_host, conn$couch_http_port, sep=",")
  url <- paste(proto, "://", conn$couch_http_host, ":", conn$couch_http_port, sep="")
  paste("CouchDB connection to: (", url, ")")
}

couch_check_status <- function(conn, expected_codes, response) {
  status_code <- response$status_code
  if(any(expected_codes == status_code)) {
    response
  } else {
    # TODO - better error handling
    simpleError("Error in response from CouchDB")
  }
}

#' @title Ping connection
#' @description Check connection to the database.
#' @param conn A couchDB connection object.
#' @export
couch_ping <- function(conn) {
  path <- couch_base_url(conn)
  expected_codes = c(200)
  result <- GET(path)
  couch_check_status(conn, expected_codes, result)
}

# Internal use
couch_store_headers_put <- function(content_type, opts) {
  c("Content-Type"=content_type);
}

# Internal use.
couch_store_headers_post <- function(content_type, opts) {
  c("Content-Type"=content_type)
}

# Couch_fetch_raw returns the entire HTTP response
# you'll need to decode the response using content()
couch_fetch_raw <- function(conn, database, key, opts=NULL) {
  path <- couch_fetch_url(conn, database, key, opts)
  expected_codes <- c(200, 300, 304)
  result <- GET(path)
  status_code <- result$status_code
  if(any(expected_codes == status_code)) {
    result
  } else {
    # TODO
    warning("Error fetching value from CouchDB");
    result
  }
}


#' @title Add attachment to document
#' @param location url of the 
#' @param revtag revision (version)
#' @param attachment file to attach
#' @param content_type Content type of the attachment (for example: "image/png")
#' @description
#' Send attachment to an existing url
#' @export
couch_attach <- function(location,revtag,attachment,content_type="image/png"){
  path <- paste0(location,"/",basename(attachment),"?rev=",revtag);
  headers <- couch_store_headers_put(content_type)
  command <- paste0('curl -vX PUT ',path,' --data-binary @',attachment, ' -H "Content-Type:',content_type,'"');
  result <- system(command);
}


#'@title Fetch a document/record.
#'@description Fetches a couch object based on the key
#'@param conn  A couchDB connection object
#'@param database  The database to connect to.
#'@param key Key of document to fetch
#'@param myOpts Additional options (not implemented in this version) 
#'@return A list object with the values from the record.
#'@export
couch_fetch <- function(conn, database, key, myOpts=NULL) {
  result <- couch_fetch_raw(conn, database, key, myOpts)
  fromJSON(content(result, as="parsed"))
}

#'@title Fetch attachment 
#'@description Gets a named attachment from a recors
#'@param conn A couchDB connection object
#'@param database The database to connecto to
#'@param key Key of document
#'@param attachment name of attachment
#'@export
couch_fetch_attachment <- function(conn,database,key,attachment){
  key <- paste0(key,"/",attachment);
  couch_fetch_raw(conn,database,key,NULL);
}

#'@title Fetch document/record from default store.
#'@description Fetches a document specified by Key from the default database on the default connection
#'@param key  The key of the document to fetch
#'@param myOpts Additional options (not implemented in this version) 
#'@return A list object with the values from the record.
#'@export
couch_fetch_default <- function(key, myOpts=NULL){
  couch_fetch(conn=couch_default_connection,database=couch_default_database,key,myOpts);  
}


#'@title New couchDB document
#'@param value  List to be converted to json for transmission or preformatted JSON string
#'@param database  The database to use
#'@param key  The key (recordname) to use for the object.
#'@description Creates a new object to to insert to the couchDB.
#'Takes either a list or a formatted json object as value
#'Any attachment to the record needs to be base64-enconded added to the list as "_attachments"
#'If key is provided this is used, null sends a key-less record to couch and the key will have to be retrieved from the reponse object.
#' @examples \dontrun{ 
#'    # This code creates a document containing a small list for storage in the "localhost" 
#'    # database with the key "testDoc".
#'    myDoc <- couch_new_object(list(a=1,b=2),"localhost","testDoc"); 
#'    
#'    #Same as above but with json entered directly (not recommended).
#'    myDoc <- couch_new_object('{"a":1,"b":2}',"localhost","testDoc")
#' } 
#'@export
couch_new_object <- function(value, database=NULL, key=NULL) {
  if (is.null(database))database <- couch_default_database;
  if(is.list(value))value<-toJSON(value);
  list(value=value, database=database, key=key, content_type="application/json")
}

#'@title Create database
#'@description Creates a new database based on the dbname.
#'@param conn a couchDB connection object
#'@param dbname the name of the database
#'@description Creates a new couchDB database on the connection provided.
#' @examples \dontrun{ 
#' #Note: this example assumes that there is a couchDB instance available on localhost
#'    myConn <- couch_http_connection("localhost");
#'    couch_create_database(myConn,"myDatabase") 
#' } 
#'@export
couch_create_database<- function(conn,dbname){
  path = paste0(couch_base_url(conn),"/",dbname);
  result <- PUT(path);
  result;
}

#'@title Store a record
#'@param conn A couchDB connection object
#'@param obj  A list formatted by calling couch_new_object
#'@param myOpts Additional options (not implemented in this version) 
#'@description Stores a record to the connection provided (database spec is in object )
#'@export
couch_store <- function(conn, obj, myOpts=NULL) { 
  path <- couch_store_url(conn, obj$database, obj$key)
  expected_codes <- c(200, 201, 204, 300);
  accept_json()
  headers <- couch_store_headers_post(obj$content_type, myOpts);
  if(is.null(obj$key)){
    result <- POST(path, body=obj$value,add_headers(headers));
  } else {
    result <- PUT(path, body=obj$value,add_headers(headers));
  }
  result
}

#'@title Store a document on the default connection.
#'@param obj A list formatted by calling couch_new_object.
#'@param myOpts Additional options (not implemented in this version) 
#'@description Stores a record on the default connection.
#'@export
couch_store_default <- function(obj, myOpts=NULL) { 
  path <- couch_store_url(couch_default_connection, obj$database, obj$key)
  expected_codes <- c(200, 201, 204, 300);
  accept_json()
  headers <- couch_store_headers_post(obj$content_type, myOpts);
  if(is.null(obj$key)){
    result <- POST(path, body=obj$value,add_headers(headers));
  } else {
    result <- PUT(path, body=obj$value,add_headers(headers));
  }
  result
}


#'@title Delete a record.
#'@description Delete a record on the connection provided
#'@param conn  A couchDB connection object
#'@param database Name of database to operate on
#'@param key key (record) to delete
#'@param myOpts Additional options (not implemented in this version) 
#'@export
couch_delete <- function(conn, database, key, myOpts=NULL) {
  path <- couch_delete_url(conn, database, key, myOpts)
  expected_codes <- c(204, 404)
  result <- DELETE(path)
  result
}


couch_mapreduce <- function(conn, query) {
  # TODO: chunked option
  # TODO: allow for storing views
  path <- couch_mapred_url(conn)
  expected_codes <- c(200)
  add_headers("ContentType: application/json")
  result <- POST(path,body=query);
  content(couch_check_status(conn, expected_codes, result))
}

#'@title List available databases.
#'@param conn A couchDB connection object
#'@description Lists the available databases on the connection provided.
#'@export
couch_list_databases <- function(conn) {
  path <- couch_list_databases_url(conn);
  expected_codes <- c(200);
  fromJSON(content(GET(path)))
}
