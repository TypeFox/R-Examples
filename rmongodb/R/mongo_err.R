#' Retrieve an connection error code from a mongo object
#' 
#' Retrieve an connection error code from a mongo object indicating the failure
#' code if mongo.create() failed.
#' 
#' @param mongo (\link{mongo}) a mongo connection object.
#' 
#' @return (integer) error code as follows:
#'  \itemize{
#' \item 0 No Error
#' \item 1 No socket - Could not create socket.
#' \item 2 Fail - An error occurred attempting to connect to socket
#' \item 3 Address fail - An error occured calling getaddrinfo().
#' \item 4 Not Master - Warning: connected to a non-master node
#' (read-only).
#' \item 5 Bad set name - given name doesn't match the replica set.
#' \item 6 No Primary - Cannot find primary in replica set - connection
#' closed.
#' \item 7 I/O error - An error occured reading or writing on the socket.
#' \item 8 Read size error - The response is not the expected length.
#' \item 9 Command failed - The command returned with 'ok' value of 0.
#' \item 10 BSON invalid - Not valid for the specified operation.
#' \item 11 BSON not finished - should not occur with R driver.
#' }
#' 
#' @seealso \code{\link{mongo.create}},\cr \link{mongo}
#' 
#' @examples
#' mongo <- mongo.create()
#' if (!mongo.is.connected(mongo)) {
#'     print("Unable to connect.  Error code:")
#'     print(mongo.get.err(mongo))
#' }
#' 
#' @export mongo.get.err
mongo.get.err <- function(mongo){
  .Call(".mongo.get.err", mongo)
}



#' Retrieve an server error code from a mongo connection object
#' 
#' Retrieve an server error record from a the MongoDB server.  This describes
#' the last error that occurs while accessing the give database. While this
#' function retrieves an error record in the form of a mongo.bson record, it
#' also sets the values returned by \code{\link{mongo.get.server.err}()} and
#' \code{\link{mongo.get.server.err.string}()}. You may find it more convenient
#' using those after calling \code{mongo.get.last.err()} rather than unpacking
#' the returned mongo.bson object.
#' 
#' 
#' @param mongo (\link{mongo}) a mongo connection object.
#' @param db (string) The name of the database for which to get the error
#' status.
#' @return NULL if no error was reported; otherwise,
#' 
#' (\link{mongo.bson}) This BSON object has the form { err : "\emph{error
#' message string}", code : \emph{error code integer} }
#' @seealso \code{\link{mongo.get.server.err}},\cr
#' \code{\link{mongo.get.server.err.string}},\cr
#' \code{\link{mongo.get.prev.err}}\cr \link{mongo}.
#' @examples
#' 
#' mongo <- mongo.create()
#' if (mongo.is.connected(mongo)) {
#' 
#'     # try adding a duplicate record when index doesn't allow this
#' 
#'     db <- "test"
#'     ns <- "test.people"
#'     mongo.index.create(mongo, ns, '{"name":1}', mongo.index.unique)
#' 
#'     buf <- mongo.bson.buffer.create()
#'     mongo.bson.buffer.append(buf, "name", "John")
#'     mongo.bson.buffer.append(buf, "age", 22L)
#'     b <- mongo.bson.from.buffer(buf)
#'     mongo.insert(mongo, ns, b);
#' 
#'     buf <- mongo.bson.buffer.create()
#'     mongo.bson.buffer.append(buf, "name", "John")
#'     mongo.bson.buffer.append(buf, "age", 27L)
#'     b <- mongo.bson.from.buffer(buf)
#'     mongo.insert(mongo, ns, b);
#' 
#'     err <- mongo.get.last.err(mongo, db)
#'     print(mongo.get.server.err(mongo))
#'     print(mongo.get.server.err.string(mongo))
#' }
#' 
#' @export mongo.get.last.err
mongo.get.last.err <- function(mongo, db){
  .Call(".mongo.get.last.err", mongo, db)
}


#' Retrieve an server error code from a mongo connection object
#' 
#' Retrieve the previous server error record from a the MongoDB server.  While
#' this function retrieves an error record in the form of a mongo.bson record,
#' it also sets the values returned by \code{\link{mongo.get.server.err}()} and
#' \code{\link{mongo.get.server.err.string}()}. You may find it more convenient
#' using those after calling \code{mongo.get.prev.err()} rather than unpacking
#' the returned mongo.bson object.
#' 
#' 
#' @param mongo (\link{mongo}) a mongo connection object.
#' @param db (string) The name of the database for which to get the error
#' status.
#' @return NULL if no error was reported; otherwise,
#' 
#' (\link{mongo.bson}) This BSON object has the form { err : "\emph{error
#' message string}", code : \emph{error code integer} }
#' @seealso \code{\link{mongo.get.server.err}},\cr
#' \code{\link{mongo.get.server.err.string}},\cr
#' \code{\link{mongo.get.last.err}}\cr \link{mongo}.
#' @examples
#' 
#' mongo <- mongo.create()
#' if (mongo.is.connected(mongo)) {
#' 
#'     # try adding a duplicate record when index doesn't allow this
#' 
#'     db <- "test"
#'     ns <- "test.people"
#'     mongo.index.create(mongo, ns, '{"name":1}', mongo.index.unique)
#' 
#'     buf <- mongo.bson.buffer.create()
#'     mongo.bson.buffer.append(buf, "name", "John")
#'     mongo.bson.buffer.append(buf, "age", 22L)
#'     b <- mongo.bson.from.buffer(buf)
#'     mongo.insert(mongo, ns, b);
#' 
#'     buf <- mongo.bson.buffer.create()
#'     mongo.bson.buffer.append(buf, "name", "John")
#'     mongo.bson.buffer.append(buf, "age", 27L)
#'     b <- mongo.bson.from.buffer(buf)
#'     mongo.insert(mongo, ns, b);
#' 
#'     # try insert again
#'     mongo.insert(mongo, ns, b);
#' 
#'     err <- mongo.get.prev.err(mongo, db)
#'     print(mongo.get.server.err(mongo))
#'     print(mongo.get.server.err.string(mongo))
#' }
#' 
#' @export mongo.get.prev.err
mongo.get.prev.err <- function(mongo, db){
  .Call(".mongo.get.prev.err", mongo, db)
}


#' Retrieve an server error code from a mongo connection object
#' 
#' Send a "reset error" command to the server, it also resets the values
#' returned by\cr \code{\link{mongo.get.server.err}()} and
#' \code{\link{mongo.get.server.err.string}()}.
#' 
#' 
#' @param mongo (\link{mongo}) a mongo connection object.
#' @param db (string) The name of the database on which to reset the error
#' status.
#' @return NULL
#' @seealso \code{\link{mongo.get.server.err}},\cr
#' \code{\link{mongo.get.server.err.string}},\cr
#' \code{\link{mongo.get.last.err}},\cr \code{\link{mongo.get.prev.err}},\cr
#' \link{mongo}.
#' @examples
#' 
#' mongo <- mongo.create()
#' if (mongo.is.connected(mongo)) {
#' 
#'     # try adding a duplicate record when index doesn't allow this
#' 
#'     db <- "test"
#'     ns <- "test.people"
#'     mongo.index.create(mongo, ns, '{"name":1}', mongo.index.unique)
#' 
#'     buf <- mongo.bson.buffer.create()
#'     mongo.bson.buffer.append(buf, "name", "John")
#'     mongo.bson.buffer.append(buf, "age", 22L)
#'     b <- mongo.bson.from.buffer(buf)
#'     mongo.insert(mongo, ns, b);
#' 
#'     buf <- mongo.bson.buffer.create()
#'     mongo.bson.buffer.append(buf, "name", "John")
#'     mongo.bson.buffer.append(buf, "age", 27L)
#'     b <- mongo.bson.from.buffer(buf)
#'     mongo.insert(mongo, ns, b);
#' 
#'     err <- mongo.get.last.err(mongo, db)
#'     print(mongo.get.server.err(mongo))
#'     print(mongo.get.server.err.string(mongo))
#'     mongo.reset.err(mongo, db)
#' }
#' 
#' @export mongo.reset.err
mongo.reset.err <- function(mongo, db){  
  .Call(".mongo.reset.err", mongo, db)
}


#' Retrieve an server error code from a mongo connection object
#' 
#' Retrieve an server error code from a mongo connection object.
#' 
#' \code{\link{mongo.find}()}, \code{\link{mongo.find.one}()},
#' \code{\link{mongo.index.create}()} set or clear this error code depending on
#' whether they are successful or not.
#' 
#' \code{\link{mongo.get.last.err}()} and \code{\link{mongo.get.prev.err}()}
#' both set or clear this error code according to what the server reports.
#' 
#' 
#' @param mongo (\link{mongo}) a mongo connection object.
#' @return (integer) Server error code
#' @seealso \code{\link{mongo.get.server.err.string}},\cr
#' \code{\link{mongo.get.last.err}},\cr \code{\link{mongo.get.prev.err}},\cr
#' \code{\link{mongo.find}},\cr \code{\link{mongo.find.one}},\cr
#' \code{\link{mongo.index.create}},\cr \link{mongo}.
#' @examples
#' 
#' mongo <- mongo.create()
#' if (mongo.is.connected(mongo)) {
#'     # construct a query containing invalid operator
#'     buf <- mongo.bson.buffer.create()
#'     mongo.bson.buffer.start.object(buf, "age")
#'     mongo.bson.buffer.append(buf, "$bad", 1L)
#'     mongo.bson.buffer.finish.object(buf)
#'     query <- mongo.bson.from.buffer(buf)
#' 
#'     result <- mongo.find.one(mongo, "test.people", query)
#'     if (is.null(result)) {
#'         print(mongo.get.server.err.string(mongo))
#'         print(mongo.get.server.err(mongo))
#'     }
#' }
#' 
#' @export mongo.get.server.err
mongo.get.server.err <- function(mongo){
  .Call(".mongo.get.server.err", mongo)
}


#' Retrieve an server error code from a mongo connection object
#' 
#' Retrieve an server error string from a mongo connection object.
#' 
#' \code{\link{mongo.find}()}, \code{\link{mongo.find.one}()},
#' \code{\link{mongo.index.create}()} set or clear this error string depending
#' on whether they are successful or not.
#' 
#' \code{\link{mongo.get.last.err}()} and \code{\link{mongo.get.prev.err}()}
#' both set or clear this error string according to what the server reports.
#' 
#' 
#' @param mongo (\link{mongo}) a mongo connection object.
#' @return (string) Server error string
#' @seealso \code{\link{mongo.get.server.err}},\cr
#' \code{\link{mongo.get.last.err}},\cr \code{\link{mongo.get.prev.err}},\cr
#' \code{\link{mongo.find}},\cr \code{\link{mongo.find.one}},\cr
#' \code{\link{mongo.index.create}},\cr \link{mongo}.
#' @examples
#' 
#' mongo <- mongo.create()
#' if (mongo.is.connected(mongo)) {
#'     # construct a query containing invalid operator
#'     buf <- mongo.bson.buffer.create()
#'     mongo.bson.buffer.start.object(buf, "age")
#'     mongo.bson.buffer.append(buf, "$bad", 1L)
#'     mongo.bson.buffer.finish.object(buf)
#'     query <- mongo.bson.from.buffer(buf)
#' 
#'     result <- mongo.find.one(mongo, "test.people", query)
#'     if (is.null(result)) {
#'         print(mongo.get.server.err(mongo))
#'         print(mongo.get.server.err.string(mongo))
#'     }
#' }
#' 
#' @export mongo.get.server.err.string
mongo.get.server.err.string <- function(mongo){
  .Call(".mongo.get.server.err.string", mongo)
}
