#' Add record to a collection
#' 
#' Add record to a collection.
#' 
#' See \url{http://www.mongodb.org/display/DOCS/Inserting}.
#' 
#' 
#' @param mongo (\link{mongo}) a mongo connection object.
#' @param ns (string) namespace of the collection to which to add the record.
#' @param b (\link{mongo.bson}) The record to add.
#' 
#' In addition, \code{b} may be a list which will be converted to a mongo.bson
#' object by \code{\link{mongo.bson.from.list}()}.
#' @return TRUE if the command was successfully sent to the server; otherwise,
#' FALSE.
#' 
#' \code{\link{mongo.get.last.err}()} may be examined to verify that the insert
#' was successful on the server if necessary.
#' @seealso \code{\link{mongo.insert.batch}},\cr \code{\link{mongo.update}},\cr
#' \code{\link{mongo.find}},\cr \code{\link{mongo.find.one}},\cr
#' \code{\link{mongo.remove}},\cr \link{mongo.bson},\cr \link{mongo}.
#' @examples
#' 
#' mongo <- mongo.create()
#' if (mongo.is.connected(mongo)) {
#'     ns <- "test.people"
#' 
#'     buf <- mongo.bson.buffer.create()
#'     mongo.bson.buffer.append(buf, "name", "Joe")
#'     mongo.bson.buffer.append(buf, "age", 22L)
#'     b <- mongo.bson.from.buffer(buf)
#'     mongo.insert(mongo, ns, b)
#' 
#'     # do the same thing in shorthand:
#'     mongo.insert(mongo, ns, list(name="Joe", age=22L))
#' }
#' 
#' @export mongo.insert
mongo.insert <- function(mongo, ns, b) {
  
  #check for mongodb connection
  if( !mongo.is.connected(mongo))
    stop("No mongoDB connection!")
  
  #validate and process input
  b <- switch( class(b),
                   "mongo.bson" = b,
                   "list" = mongo.bson.from.list(b),
                   "character" = mongo.bson.from.JSON(b) )

  .Call(".mongo.insert", mongo, ns, b)
}



#' Add multiple records to a collection
#' 
#' Add multiple records to a collection.  This function eliminates some network
#' traffic and server overhead by sending all the records in a single message.
#' 
#' See \url{http://www.mongodb.org/display/DOCS/Inserting}.
#' 
#' 
#' @param mongo (\link{mongo}) a mongo connection object.
#' @param ns (string) namespace of the collection to which to add the record.
#' @param lst A list of (\link{mongo.bson}) records to add.
#' @return TRUE if the command was successfully sent to the server; otherwise,
#' FALSE.
#' 
#' \code{\link{mongo.get.last.err}()} may be examined to verify that the insert
#' was successful on the server if necessary.
#' @seealso \code{\link{mongo.insert}},\cr \code{\link{mongo.update}},\cr
#' \code{\link{mongo.find}},\cr \code{\link{mongo.find.one}},\cr
#' \code{\link{mongo.remove}},\cr \link{mongo.bson},\cr \link{mongo}.
#' @examples
#' 
#' mongo <- mongo.create()
#' if (mongo.is.connected(mongo)) {
#'     ns <- "test.people"
#' 
#'     buf <- mongo.bson.buffer.create()
#'     mongo.bson.buffer.append(buf, "name", "Dave")
#'     mongo.bson.buffer.append(buf, "age", 27L)
#'     x <- mongo.bson.from.buffer(buf)
#' 
#'     buf <- mongo.bson.buffer.create()
#'     mongo.bson.buffer.append(buf, "name", "Fred")
#'     mongo.bson.buffer.append(buf, "age", 31L)
#'     y <- mongo.bson.from.buffer(buf)
#' 
#'     buf <- mongo.bson.buffer.create()
#'     mongo.bson.buffer.append(buf, "name", "Silvia")
#'     mongo.bson.buffer.append(buf, "city", 24L)
#'     z <- mongo.bson.from.buffer(buf)
#'     mongo.insert.batch(mongo, ns, list(x, y, z))
#' }
#' 
#' @export mongo.insert.batch
mongo.insert.batch <- function(mongo, ns, lst){
  
  #check for mongodb connection
  if( !mongo.is.connected(mongo))
    stop("No mongoDB connection!")

  res <- .Call(".mongo.insert.batch", mongo, ns, lst)
  
  if( res == FALSE ){
    warning(mongo.get.server.err.string(mongo))
  }
  
  return(res)
}



#' mongo.update() flag constant for an upsert
#' 
#' Flag to \code{\link{mongo.update}()} (1L): insert ObjNew into the database
#' if no record matching criteria is found.
#' 
#' 
#' @return 1L
#' @seealso \code{\link{mongo.update}},\cr \code{\link{mongo.update.multi}},\cr
#' \code{\link{mongo.update.basic}}.
#' @export mongo.update.upsert
mongo.update.upsert <- 1L


#' mongo.update() flag constant for updating multiple records
#' 
#' Flag to \code{\link{mongo.update}()} (2L): Update multiple records rather
#' than just the first one matched by criteria.
#' 
#' 
#' @return 2L
#' @seealso \code{\link{mongo.update}},\cr
#' \code{\link{mongo.update.upsert}},\cr \code{\link{mongo.update.basic}}.
#' @export mongo.update.multi
mongo.update.multi  <- 2L


#' mongo.update() flag constant for performing a basic update
#' 
#' Flag to \code{\link{mongo.update}()} (4L): Perform a basic update.
#' 
#' 
#' @return 4L
#' @seealso \code{\link{mongo.update}},\cr \code{\link{mongo.update.multi}}\cr
#' \code{\link{mongo.update.upsert}}
#' @export mongo.update.basic
mongo.update.basic  <- 4L



#' Perform an update on a collection
#' 
#' Perform an update on a collection.
#' 
#' See \url{http://www.mongodb.org/display/DOCS/Updating}.
#' 
#' 
#' @param mongo (\link{mongo}) a mongo connection object.
#' @param ns (string) namespace of the collection to which to update.
#' @param criteria (\link{mongo.bson}) The criteria with which to match records
#' that are to be updated.
#' 
#' Alternately, \code{criteria} may be a list which will be converted to a
#' mongo.bson object by \code{\link{mongo.bson.from.list}()}.
#' 
#' Alternately, \code{criteria} may be a valid JSON character string which will be converted to a
#' mongo.bson object by \code{\link{mongo.bson.from.JSON}()}.
#' @param objNew (\link{mongo.bson}) The replacement object.
#' 
#' Alternately, \code{objNew} may be a list which will be converted to a
#' mongo.bson object by \code{\link{mongo.bson.from.list}()}.
#' 
#' Alternately, \code{objNew} may be a valid JSON character string which will be converted to a
#' mongo.bson object by \code{\link{mongo.bson.from.JSON}()}.
#' @param flags (integer vector) A list of optional flags governing the
#' operation: \itemize{ \item\code{\link{mongo.update.upsert}}: insert ObjNew
#' into the database if no record matching criteria is found.
#' \item\code{\link{mongo.update.multi}}: update multiple records rather than
#' just the first one matched by criteria.
#' \item\code{\link{mongo.update.basic}}: Perform a basic update.  }
#' @seealso \link{mongo},\cr \link{mongo.bson},\cr
#' \code{\link{mongo.insert}},\cr \code{\link{mongo.find}},\cr
#' \code{\link{mongo.find.one}},\cr \code{\link{mongo.remove}}.
#' @examples
#' 
#' mongo <- mongo.create()
#' if (mongo.is.connected(mongo)) {
#'     ns <- "test.people"
#' 
#'     buf <- mongo.bson.buffer.create()
#'     mongo.bson.buffer.append(buf, "name", "Joe")
#'     criteria <- mongo.bson.from.buffer(buf)
#' 
#'     buf <- mongo.bson.buffer.create()
#'     mongo.bson.buffer.start.object(buf, "$inc")
#'     mongo.bson.buffer.append(buf, "age", 1L)
#'     mongo.bson.buffer.finish.object(buf)
#'     objNew <- mongo.bson.from.buffer(buf)
#' 
#'     # increment the age field of the first record matching name "Joe"
#'     mongo.update(mongo, ns, criteria, objNew)
#' 
#'     buf <- mongo.bson.buffer.create()
#'     mongo.bson.buffer.append(buf, "name", "Jeff")
#'     criteria <- mongo.bson.from.buffer(buf)
#' 
#'     buf <- mongo.bson.buffer.create()
#'     mongo.bson.buffer.append(buf, "name", "Jeff")
#'     mongo.bson.buffer.append(buf, "age", 27L)
#'     objNew <- mongo.bson.from.buffer(buf)
#' 
#'     # update the entire record to { name: "Jeff", age: 27 }
#'     # where name equals "Jeff"
#'     # if such a record exists; otherwise, insert this as a new reord
#'     mongo.update(mongo, ns, criteria, objNew,
#'         mongo.update.upsert)
#' 
#'     # do a shorthand update:
#'     mongo.update(mongo, ns, list(name="John"), list(name="John", age=25))
#' }
#' 
#' @export mongo.update
mongo.update <- function(mongo, ns, criteria, objNew, flags=0L) {
  
  #check for mongodb connection
  if( !mongo.is.connected(mongo))
    stop("No mongoDB connection!")
  
  #validate and process input
  criteria <- switch( class(criteria),
                   "mongo.bson" = criteria,
                   "list" = mongo.bson.from.list(criteria),
                   "character" = mongo.bson.from.JSON(criteria) )
  objNew <- switch( class(objNew),
                  "mongo.bson" = objNew,
                  "list" = mongo.bson.from.list(objNew),
                  "character" = mongo.bson.from.JSON(objNew) )

  .Call(".mongo.update", mongo, ns, criteria, objNew, flags)
}



#' Remove records from a collection
#' 
#' Remove all records from a collection that match a given criteria.
#' 
#' See \url{http://www.mongodb.org/display/DOCS/Removing}.
#' 
#' 
#' @param mongo (\link{mongo}) a mongo connection object.
#' @param ns (string) namespace of the collection from which to remove records.
#' @param criteria (\link{mongo.bson}) The criteria with which to match records
#' that are to be removed. The default of mongo.bson.empty() will cause
#' \emph{all} records in the given collection to be removed.
#' 
#' Alternately, \code{criteria} may be a list which will be converted to a
#' mongo.bson object by \code{\link{mongo.bson.from.list}()}.
#' 
#' Alternately, \code{criteria} may be a valid JSON character string which will be converted to a
#' mongo.bson object by \code{\link{mongo.bson.from.JSON}()}.
#' @seealso \link{mongo},\cr \link{mongo.bson},\cr
#' \code{\link{mongo.insert}},\cr \code{\link{mongo.update}},\cr
#' \code{\link{mongo.find}},\cr \code{\link{mongo.find.one}}.
#' @examples
#' 
#' mongo <- mongo.create()
#' if (mongo.is.connected(mongo)) {
#'     buf <- mongo.bson.buffer.create()
#'     mongo.bson.buffer.append(buf, "name", "Jeff")
#'     criteria <- mongo.bson.from.buffer(buf)
#' 
#'     # remove all records where name is "Jeff"
#'     # from collection people in database test
#'     mongo.remove(mongo, "test.people", criteria)
#' 
#'     # remove all records from collection cars in database test
#'     mongo.remove(mongo, "test.cars")
#' 
#'     # shorthand: remove all records where name is "Fred"
#'     mongo.remove(mongo, "test.people", list(name="Fred"))
#' }
#' 
#' @export mongo.remove
mongo.remove <- function(mongo, ns, criteria=mongo.bson.empty()) {
  
  #check for mongodb connection
  if( !mongo.is.connected(mongo))
    stop("No mongoDB connection!")
  
  #validate and process input
  criteria <- switch( class(criteria),
                      "mongo.bson" = criteria,
                      "list" = mongo.bson.from.list(criteria),
                      "character" = mongo.bson.from.JSON(criteria) )
  
  .Call(".mongo.remove", mongo, ns, criteria)
}

