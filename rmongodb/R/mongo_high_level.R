#' Get a vector of distinct values for keys in a collection
#'
#' Get a vector of distinct values for keys in a collection.
#'
#' See
#' \url{http://www.mongodb.org/display/DOCS/Aggregation#Aggregation-Distinct}.
#'
#'
#' @param mongo (\link{mongo}) A mongo connection object.
#' @param ns (string) The namespace of the collection in which to find distinct
#' keys.
#' @param key (string) The name of the key field for which to get distinct
#' values.
#' @param query \link{mongo.bson} An optional query to restrict the returned
#' values.
#'
#' @return vector of distinct values or NULL if the command failed.
#'
#' (vector) The result set of distinct keys.
#' @seealso \code{\link{mongo.command}},\cr
#' \code{\link{mongo.simple.command}},\cr \code{\link{mongo.find}},\cr
#' \link{mongo}.
#'
#' @examples
#' mongo <- mongo.create()
#' if (mongo.is.connected(mongo)) {
#'     keys <- mongo.distinct(mongo, "test.people", "name")
#'     print(keys)
#' }
#'
#' @aliases mongo.get.values
#' @export mongo.get.values
#' @export mongo.distinct
mongo.distinct <- function(mongo, ns, key, query=mongo.bson.empty()) {

  ns_parsed <- mongo.parse.ns(ns)
  db <- ns_parsed$db
  collection <- ns_parsed$collection
  if( is.null(db) || is.null(collection) ){
    stop("Wrong namespace (ns).")
  }

  b <- mongo.command(mongo, db, list(distinct=collection, key=key, query=query))

  if (!is.null(b)){
    b <- mongo.bson.value(b, "values")
    if(length(b)==0)
      warning("No values - probably wrong key!")
    return(b)
  } else{
    warning( mongo.get.server.err.string(mongo) )
    return(NULL)
  }
}

mongo.get.values <- mongo.distinct





#' Aggregation pipeline
#'
#' Aggregation pipeline
#'
#' See
#' \url{http://docs.mongodb.org/manual/reference/command/aggregate/}
#' \url{http://docs.mongodb.org/manual/core/aggregation-pipeline/}.
#'
#' @param mongo (\link{mongo}) A mongo connection object.
#' @param ns (string) The namespace of the collection in which to find distinct
#' keys.
#' @param pipeline (\link{list} of \link{mongo.bson} objects) representing aggregation query pipeline.
#' Alternately, \code{pipeline} may be a \link{list} of \link{list} which will be converted to a mongo.bson list object by
#' \code{\link{mongo.bson.from.list}()}.
#'
#' Alternately, \code{pipeline} may be a \link{list} of valid JSON \link{character} strings which will be converted to a
#' mongo.bson object by \code{\link{mongo.bson.from.JSON}()}.
#' @param explain (\link{logical}) Optional, MongoDB 2.6+. Specifies to return the information on the processing of the pipeline. References above.
#' @param allowDiskUse (\link{logical}) Optional, MongoDB 2.6+. Enables writing to temporary files. When set to true, aggregation stages can write data to the _tmp subdirectory in the dbPath directory.
#' @param cursor (\link{mongo.bson}) Optional, MongoDB 2.6+. Specify a document that contains options that control the creation of the cursor object.
#' @param ... Arguments to be passed to methods, such as \link{mongo.bson.to.list}, \link{fromJSON}
#' Unfortunately, current underlying mongo-c-driver can return BSON from aggreagation camand. Cursors are not supported.
#'
#' Alternately, \code{cursor} may be a list which will be converted to a
#' mongo.bson object by \code{\link{mongo.bson.from.list}()}.
#'
#' Alternately, \code{cursor} may be a valid JSON character string which will be converted to mongo.bson object by \code{\link{mongo.bson.from.JSON}()}.
#'
#' @return NULL if the command failed.  \code{\link{mongo.get.err}()} may be
#' MONGO_COMMAND_FAILED.
#'
#' \link{mongo.bson} The result of aggregation.
#' @seealso \code{\link{mongo.command}},\cr
#' \code{\link{mongo.simple.command}},\cr \code{\link{mongo.find}},\cr
#' \link{mongo}.
#'
#' @examples
#' # using the zips example data set
#' mongo <- mongo.create()
#' # insert some example data
#' data(zips)
#' colnames(zips)[5] <- "orig_id"
#' ziplist <- list()
#' ziplist <- apply( zips, 1, function(x) c( ziplist, x ) )
#' res <- lapply( ziplist, function(x) mongo.bson.from.list(x) )
#' if (mongo.is.connected(mongo)) {
#'     mongo.insert.batch(mongo, "test.zips", res )
#'     pipe_1 <- mongo.bson.from.JSON('{"$group":{"_id":"$state", "totalPop":{"$sum":"$pop"}}}')
#'     cmd_list <- list(pipe_1)
#'     res <- mongo.aggregation(mongo, "test.zips", cmd_list)
#' }
#' mongo.destroy(mongo)
#'
#' @export mongo.aggregation
mongo.aggregation <- function(mongo, ns, pipeline, explain = NULL, allowDiskUse = NULL, cursor = NULL, ...)
{
  if(!inherits(pipeline, what = 'list')) stop("pipeline should be a list!")
  ns_parsed <- mongo.parse.ns(ns)
  db <- ns_parsed$db
  collection <- ns_parsed$collection
  if( is.null(db) || is.null(collection) ){
    stop("Wrong namespace (ns).")
  }
  command <- list()
  command[['aggregate']] <- collection
  command[['pipeline']] <- lapply(pipeline, mongo.list.from.argument, ...)
  # New parameters for MongoDB 2.6+
  if(is.logical(explain)) command[['explain']] <- explain
  if(is.logical(allowDiskUse)) command[['allowDiskUse']] <- allowDiskUse
  if(!is.null(cursor)) command[['cursor']] <- mongo.list.from.argument(cursor)
  # Make final bson command
  command <- mongo.bson.from.list(command)

  res <- mongo.command(mongo, db, command)

  if( is.null(res) ){
    stop(paste("mongoDB error: ", mongo.get.err(mongo), ". Please check ?mongo.get.err for more details.", sep=""))
  }

  return(res)
}
