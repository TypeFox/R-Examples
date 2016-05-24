#' Find one record in a collection
#'
#' Find the first record in a collection that matches a given query.
#'
#' This is a simplified version of mongo.find() which eliminates the need to
#' step through returned records with a cursor.
#'
#' See \url{http://www.mongodb.org/display/DOCS/Querying}.
#'
#'
#' @param mongo (\link{mongo}) A mongo connection object.
#' @param ns (string) The namespace of the collection from in which to find a
#' record.
#' @param query (\link{mongo.bson}) The criteria with which to match the record
#' that is to be found. The default of mongo.bson.empty() will cause the the
#' very first record in the collection to be returned.
#'
#' Alternately, \code{query} may be a list which will be converted to a
#' mongo.bson object by \code{\link{mongo.bson.from.list}()}.
#'
#' Alternately, \code{query} may be a valid JSON character string which will be converted to a
#' mongo.bson object by \code{\link{mongo.bson.from.JSON}()}.
#' @param fields (\link{mongo.bson}) The desired fields which are to be
#' returned frtom the matching record.  The default of mongo.bson.empty() will
#' cause all fields of the matching record to be returned; however, specific
#' fields may be specified to cut down on network traffic and memory overhead.
#'
#' Alternately, \code{fields} may be a list which will be converted to a
#' mongo.bson object by \code{\link{mongo.bson.from.list}()}.
#'
#' Alternately, \code{fields} may be a valid JSON character string which will be converted to a
#' mongo.bson object by \code{\link{mongo.bson.from.JSON}()}.
#'
#' @return NULL if no record matching the criteria is found; otherwise,
#'
#' (\link{mongo.bson}) The matching record/fields.
#'
#' Note that NULL may also be returned if a database error occurred (when a
#' badly formed query is used, for example). \code{\link{mongo.get.server.err}}
#' and \code{\link{mongo.get.server.err.string}} may be examined in that case.
#'
#' @seealso \code{\link{mongo.find}},\cr \code{\link{mongo.index.create}},\cr
#' \code{\link{mongo.insert}},\cr \code{\link{mongo.update}},\cr
#' \code{\link{mongo.remove}},\cr \link{mongo},\cr \link{mongo.bson}.
#'
#' @examples
#' mongo <- mongo.create()
#' if (mongo.is.connected(mongo)) {
#'     buf <- mongo.bson.buffer.create()
#'     mongo.bson.buffer.append(buf, "name", "Jeff")
#'     query <- mongo.bson.from.buffer(buf)
#'
#'     # find the first record where name is "Jeff"\
#'     #    in collection people of database test
#'     b <- mongo.find.one(mongo, "test.people", query)
#'     if (!is.null(b))
#'         print(b)
#'
#'     buf <- mongo.bson.buffer.create()
#'     mongo.bson.buffer.append(buf, "_id", 1L)
#'     mongo.bson.buffer.append(buf, "age", 1L)
#'     fields <- mongo.bson.from.buffer(buf)
#'
#'     # find the first record where name is "Jeff"
#'     #    in collection people of database test
#'     # return only the _id and age fields of the matched record
#'     b <- mongo.find.one(mongo, "test.people", query, fields)
#'     if (!is.null(b))
#'         print(b)
#'
#'     # find the first record in collection cars of database test
#'     have.car <- !is.null(mongo.find.one(mongo, "test.cars"))
#'
#'     # shorthand using a list:
#'     b <- mongo.find.one(mongo, "test.people", list(name="Jose"))
#' }
#'
#' @export mongo.find.one
#' @export mongo.findOne
#'
#' @aliases mongo.findOne
mongo.find.one <- function(mongo, ns, query=mongo.bson.empty(), fields=mongo.bson.empty()) {

  #check for mongodb connection
  if( !mongo.is.connected(mongo))
    stop("No mongoDB connection!")

  #validate and process input
  query <- mongo.bson.from.argument(query)
  fields <- mongo.bson.from.argument(fields)
  .Call(".mongo.find.one", mongo, ns, query, fields)
}

mongo.findOne <- mongo.find.one



#' mongo.find flag constant - cursor tailable
#'
#' \code{\link{mongo.find}()} flag constant - cursor tailable.
#'
#'
#' @return 2L
#' @export mongo.find.cursor.tailable
mongo.find.cursor.tailable   <- 2L


#' mongo.find flag constant - slave ok
#'
#' \code{\link{mongo.find}()} flag constant - slave ok.
#'
#'
#' @return 4L
#' @export mongo.find.slave.ok
mongo.find.slave.ok          <- 4L


#' mongo.find flag constant - oplog replay
#'
#' \code{\link{mongo.find}()} flag constant - oplog replay.
#'
#'
#' @return 8L
#' @export mongo.find.oplog.replay
mongo.find.oplog.replay      <- 8L


#' mongo.find flag constant - no cursor timeout
#'
#' \code{\link{mongo.find}()} flag constant - no cursor timeout.
#'
#'
#' @return 16L
#' @export mongo.find.no.cursor.timeout
mongo.find.no.cursor.timeout <- 16L


#' mongo.find flag constant - await data
#'
#' \code{\link{mongo.find}()} flag constant - await data.
#'
#'
#' @return 32L
#' @export mongo.find.await.data
mongo.find.await.data        <- 32L


#' mongo.find flag constant - exhaust
#'
#' \code{\link{mongo.find}()} flag constant - exhaust.
#'
#'
#' @return 64L
#' @export mongo.find.exhaust
mongo.find.exhaust           <- 64L


#' mongo.find flag constant - partial results
#'
#' \code{\link{mongo.find}()} flag constant - partial results.
#'
#'
#' @return 128L
#' @export mongo.find.partial.results
mongo.find.partial.results   <- 128L



#' Find records in a collection
#'
#' Find records in a collection that match a given query.
#'
#' See \url{http://www.mongodb.org/display/DOCS/Querying}.
#'
#'
#' @param mongo (\link{mongo}) a mongo connection object.
#' @param ns (string) namespace of the collection from which to find records.
#' @param query (\link{mongo.bson}) The criteria with which to match the
#' records to be found.  The default of mongo.bson.empty() will cause the the
#' very first record in the collection to be returned.
#'
#' Alternately, \code{query} may be a list which will be converted to a
#' mongo.bson object by \code{\link{mongo.bson.from.list}()}.
#'
#' Alternately, \code{query} may be a valid JSON character string which will be converted to a
#' mongo.bson object by \code{\link{mongo.bson.from.JSON}()}.
#' @param sort (\link{mongo.bson}) The desired fields by which to sort the
#' returned records. The default of mongo.bson.empty() indicates that no
#' special sorting is to be done; the records will come back in the order that
#' indexes locate them.
#'
#' Alternately, \code{sort} may be a list which will be converted to a
#' mongo.bson object by \code{\link{mongo.bson.from.list}()}.
#'
#' Alternately, \code{sort} may be a valid JSON character string which will be converted to a
#' mongo.bson object by \code{\link{mongo.bson.from.JSON}()}.
#' @param fields (\link{mongo.bson}) The desired fields which are to be
#' returned from the matching record.  The default of mongo.bson.empty() will
#' cause all fields of the matching record to be returned; however, specific
#' fields may be specified to cut down on network traffic and memory overhead.
#'
#' Alternately, \code{fields} may be a list which will be converted to a
#' mongo.bson object by \code{\link{mongo.bson.from.list}()}.
#'
#' Alternately, \code{fields} may be a valid JSON character string which will be converted to a
#' mongo.bson object by \code{\link{mongo.bson.from.JSON}()}.
#' @param limit (as.integer) The maximum number of records to be returned. A
#' limit of 0L will return all matching records not skipped.
#' @param skip (as.integer) The number of matching records to skip before
#' returning subsequent matching records.
#' @param options (integer vector) Flags governing the requested operation as
#' follows: \itemize{ \item\link{mongo.find.cursor.tailable}
#' \item\link{mongo.find.slave.ok} \item\link{mongo.find.oplog.replay}
#' \item\link{mongo.find.no.cursor.timeout} \item\link{mongo.find.await.data}
#' \item\link{mongo.find.exhaust} \item\link{mongo.find.partial.results} }
#' @return (\link{mongo.cursor}) An object of class "mongo.cursor" which is
#' used to step through the matching records.
#'
#' Note that an empty cursor will be returned if a database error occurred.\cr
#' \code{\link{mongo.get.server.err}()} and
#' \code{\link{mongo.get.server.err.string}()} may be examined in that case.
#' @seealso \code{\link{mongo.cursor}},\cr \code{\link{mongo.cursor.next}},\cr
#' \code{\link{mongo.cursor.value}},\cr \code{\link{mongo.find.one}},\cr
#' \code{\link{mongo.insert}},\cr \code{\link{mongo.index.create}},\cr
#' \code{\link{mongo.update}},\cr \code{\link{mongo.remove}},\cr
#' \link{mongo},\cr \link{mongo.bson}.
#' @examples
#'
#' mongo <- mongo.create()
#' if (mongo.is.connected(mongo)) {
#'     buf <- mongo.bson.buffer.create()
#'     mongo.bson.buffer.append(buf, "age", 18L)
#'     query <- mongo.bson.from.buffer(buf)
#'
#'     # Find the first 100 records
#'     #    in collection people of database test where age == 18
#'     cursor <- mongo.find(mongo, "test.people", query, limit=100L)
#'     # Step though the matching records and display them
#'     while (mongo.cursor.next(cursor))
#'         print(mongo.cursor.value(cursor))
#'     mongo.cursor.destroy(cursor)
#'
#'
#'     # shorthand: find all records where age=32, sorted by name,
#'     # and only return the name & address fields:
#'     cursor <- mongo.find(mongo, "test.people", list(age=32),
#'                          list(name=1L), list(name=1L, address=1L))
#' }
#'
#' @export mongo.find
mongo.find <- function(mongo, ns,
                       query = mongo.bson.empty(),
                       sort = mongo.bson.empty(),
                       fields = mongo.bson.empty(),
                       limit = 0L, skip = 0L, options = 0L) {

  #check for mongodb connection
  if( !mongo.is.connected(mongo))
    stop("No mongoDB connection!")

  #validate and process input
  query <- mongo.bson.from.argument(query)
  sort <- mongo.bson.from.argument(sort)
  fields <- mongo.bson.from.argument(fields)
  .Call(".mongo.find", mongo, ns, query, sort, fields, limit, skip, options)
}



#' Find records in a collection and returns one R data frame object
#'
#' Find records in a collection that match a given query and return an R data
#' frame object.
#'
#' See \url{http://www.mongodb.org/display/DOCS/Querying}.
#'
#'
#' @param mongo (\link{mongo}) a mongo connection object.
#' @param ns (string) namespace of the collection from which to find records.
#' @param query (\link{mongo.bson}) The criteria with which to match the
#' records to be found.  The default of mongo.bson.empty() will cause the the
#' very first record in the collection to be returned.
#'
#' Alternately, \code{query} may be a list which will be converted to a
#' mongo.bson object by \code{\link{mongo.bson.from.list}()}.
#'
#' Alternately, \code{query} may be a valid JSON character string which will be converted to a
#' mongo.bson object by \code{\link{mongo.bson.from.JSON}()}.
#' @param sort (\link{mongo.bson}) The desired fields by which to sort the
#' returned records. The default of mongo.bson.empty() indicates that no
#' special sorting is to be done; the records will come back in the order that
#' indexes locate them.
#'
#' Alternately, \code{sort} may be a list which will be converted to a
#' mongo.bson object by \code{\link{mongo.bson.from.list}()}.
#'
#' Alternately, \code{sort} may be a valid JSON character string which will be converted to a
#' mongo.bson object by \code{\link{mongo.bson.from.JSON}()}.
#' @param fields (\link{mongo.bson}) The desired fields which are to be
#' returned from the matching record.  The default of mongo.bson.empty() will
#' cause all fields of the matching record to be returned; however, specific
#' fields may be specified to cut down on network traffic and memory overhead.
#'
#' Alternately, \code{fields} may be a list which will be converted to a
#' mongo.bson object by \code{\link{mongo.bson.from.list}()}.
#'
#' Alternately, \code{fields} may be a valid JSON character string which will be converted to a
#' mongo.bson object by \code{\link{mongo.bson.from.JSON}()}.
#' @param limit (as.integer) The maximum number of records to be returned. A
#' limit of 0L will return all matching records not skipped.
#' @param skip (as.integer) The number of matching records to skip before
#' returning subsequent matching records.
#' @param options (integer vector) Flags governing the requested operation as
#' follows: \itemize{ \item\link{mongo.find.cursor.tailable}
#' \item\link{mongo.find.slave.ok} \item\link{mongo.find.oplog.replay}
#' \item\link{mongo.find.no.cursor.timeout} \item\link{mongo.find.await.data}
#' \item\link{mongo.find.exhaust} \item\link{mongo.find.partial.results} }
#' @param data.frame (boolean) If TRUE the result will be an \link{data.frame}
#' object, if FALSE it will be an \link{list} object. Due to NoSQL in
#' mongodb in most cases a data.frame object will not work!
#' @param mongo.oid2character (boolean) If TRUE monogo_oids will be converted to characters.
#' @param ...  optional arguments to \link{as.data.frame}
#'
#' @return An R data frame object.
#'
#' @seealso \code{\link{mongo.find.one}},\cr \code{\link{mongo.insert}},\cr
#' \code{\link{mongo.index.create}},\cr \code{\link{mongo.update}},\cr
#' \code{\link{mongo.remove}},\cr \link{mongo}.
#'
#' @examples
#' mongo <- mongo.create()
#' if (mongo.is.connected(mongo)) {
#'     buf <- mongo.bson.buffer.create()
#'     mongo.bson.buffer.append(buf, "age", 22L)
#'     query <- mongo.bson.from.buffer(buf)
#'
#'     # Find the first 100 records
#'     #    in collection people of database test where age == 22
#'     mongo.find.all(mongo, "test.people", query, limit=100L)
#'
#'
#'     # shorthand: find all records where age=22, sorted by name,
#'     # and only return the name & address fields:
#'     mongo.find.all(mongo, "test.people", list(age=22),
#'                          list(name=1L), list(name=1L, address=1L))
#' }
#'
#' @aliases mongo.find.batch
#'
#' @export mongo.find.batch
#' @export mongo.find.all
mongo.find.all <- function(mongo, ns,
                           query = mongo.bson.empty(), sort = mongo.bson.empty(),
                           fields = mongo.bson.empty(), limit = 0L, skip = 0L,
                           options = 0L,
                           data.frame = FALSE, mongo.oid2character = TRUE, ...) {

  #check for mongodb connection
  if( !mongo.is.connected(mongo))
    stop("No mongoDB connection!")

  if(data.frame==T & mongo.oid2character == F)
    warning("You won't get correct id in your data.frame if you don't set mongo.oid2character to TRUE")

  cursor <- mongo.find(mongo, ns, query = query, sort = sort, fields = fields, limit = limit, skip = skip, options = options)

  # Step though the matching records
  temp <- mongo.cursor.to.list(cursor, ...)

  if(mongo.oid2character){
    temp <- lapply(temp, function(x) {
      lapply(x, function(y){
        if(inherits(y, "mongo.oid"))
          as.character.mongo.oid(y)
        else
          y
      }
      )
    }
    )
  }


  if(data.frame){

    temp <- do.call(rbind, lapply(temp, function(x) data.frame(x, stringsAsFactors = F)) )

    if(!is.null(temp)) {
      colnames <- colnames(temp)
      select <- (colnames == "X_id")
      colnames[select] <- "_id"
      colnames(temp) <- colnames
    }
  }



  return(temp)

}
mongo.find.batch <- mongo.find.all




#' Count records in a collection
#'
#' Count the number of records in a collection that match a query See
#' \url{http://www.mongodb.org/display/DOCS/Indexes}.
#'
#'
#' @param mongo (\link{mongo}) A mongo connection object.
#' @param ns (string) The namespace of the collection in which to add count
#' records.
#' @param query \link{mongo.bson} The criteria with which to match records that
#' are to be counted.  The default of mongo.bson.empty() matches all records in
#' the collection.
#'
#' Alternately, \code{query} may be a list which will be converted to a
#' mongo.bson object by \code{\link{mongo.bson.from.list}()}.
#'
#' Alternately, \code{query} may be a valid JSON character string which will be converted to a
#' mongo.bson object by \code{\link{mongo.bson.from.JSON}()}.
#' @return (double) The number of matching records.
#' @seealso \code{\link{mongo.find}},\cr \code{\link{mongo.find.one}},\cr
#' \code{\link{mongo.insert}},\cr \code{\link{mongo.update}},\cr
#' \code{\link{mongo.remove}},\cr \link{mongo},\cr \link{mongo.bson}.
#' @examples
#'
#' mongo <- mongo.create()
#' if (mongo.is.connected(mongo)) {
#'     # Count the number of records in collection people of database test
#'     people.count <- mongo.count(mongo, "test.people")
#'     print("total people")
#'     print(people.count)
#'
#'     buf <- mongo.bson.buffer.create()
#'     mongo.bson.buffer.append(buf, "age", 21L)
#'     query <- mongo.bson.from.buffer(buf)
#'
#'     # Count the number of records in collection people of database test
#'     # where age == 21
#'     just.legal.count <- mongo.count(mongo, "test.people", query)
#'     print("people of age 21")
#'     print(just.legal.count)
#'
#'     buf <- mongo.bson.buffer.create()
#'     mongo.bson.buffer.start.object(buf, "age")
#'     mongo.bson.buffer.append(buf, "$gte", 21L)
#'     mongo.bson.buffer.finish.object(buf)
#'     query <- mongo.bson.from.buffer(buf)
#'
#'     # Count the number of records in collection people of database test
#'     # where age >= 21
#'     total.legal.count <- mongo.count(mongo, "test.people", query)
#'     print("people of age 21 or greater")
#'     print(total.legal.count)
#'
#'     # shorthand using a list:
#'     ford.count <- mongo.count(mongo, "test.cars", list(make="Ford"))
#' }
#'
#' @export mongo.count
mongo.count <- function(mongo, ns, query=mongo.bson.empty()) {

  #check for mongodb connection
  if( !mongo.is.connected(mongo))
    stop("No mongoDB connection!")

  #validate and process input
  query <- mongo.bson.from.argument(query)
  .Call(".mongo.count", mongo, ns, query)
}

