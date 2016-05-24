

#' The mongo.bson.buffer class
#' 
#' Objects of class "mongo.bson.buffer" are used to build BSON documents
#' (\link{mongo.bson} objects).
#' 
#' There are many functions for appending data into a mongo.bson.buffer
#' object.\cr See \code{\link{mongo.bson.buffer.append}()} for a list of those
#' functions.
#' 
#' After constructing your object in the buffer,
#' \code{\link{mongo.bson.from.buffer}()} may be used to turn the buffer into a
#' mongo.bson object.
#' 
#' mongo.bson.buffer objects have "mongo.bson.buffer" as their class and
#' contain an externally managed pointer to the actual document data buffer.
#' This pointer is stored in the "mongo.bson.buffer" attribute of the object.
#' 
#' 
#' @name mongo.bson.buffer
#' @docType class
#' @seealso \link{mongo.bson},\cr \code{\link{mongo.bson.buffer.size}},\cr
#' \code{\link{mongo.bson.from.buffer}},\cr
#' \code{\link{mongo.bson.buffer.append}},\cr
#' \code{\link{mongo.bson.buffer.start.object}},\cr
#' \code{\link{mongo.bson.buffer.start.array}},\cr
#' \code{\link{mongo.bson.buffer.finish.object}}.
#' @examples
#' 
#' buf <- mongo.bson.buffer.create()
#' mongo.bson.buffer.append(buf, "make", "Ford")
#' mongo.bson.buffer.append(buf, "model", "Mustang")
#' mongo.bson.buffer.append.int(buf, "year", 1968)
#' b <- mongo.bson.from.buffer(buf)
#' 
NULL





#' The mongo.bson.iterator class
#' 
#' Objects of class "mongo.bson.iterator" are used to iterate through BSON
#' documents as stored in \link{mongo.bson} objects.
#' 
#' mongo.bson.iterator objects have "mongo.bson.iterator" as their class and
#' contain an externally managed pointer to the actual document data. This
#' pointer is stored in the "mongo.bson.iterator" attribute of the object.
#' 
#' 
#' @name mongo.bson.iterator
#' @docType class
#' @seealso \code{\link{mongo.bson.iterator.create}},\cr
#' \code{\link{mongo.bson.find}},\cr \code{\link{mongo.bson.iterator.next}},\cr
#' \code{\link{mongo.bson.iterator.key}},\cr
#' \code{\link{mongo.bson.iterator.value}},\cr \link{mongo.bson}.
#' @examples
#' 
#' b <- mongo.bson.from.list(list(name="Joy", age=35, city="Ontario"))
#' # b is of class "mongo.bson"
#' iter <- mongo.bson.iterator.create(b)
#' while (mongo.bson.iterator.next(iter))
#'     print(mongo.bson.iterator.value(iter))
#' 
NULL





#' The mongo.bson class
#' 
#' Objects of class "mongo.bson" are used to store BSON documents. BSON is the
#' form that MongoDB uses to store documents in its database. MongoDB network
#' traffic also uses BSON in messages.
#' 
#' See \url{http://www.mongodb.org/display/DOCS/BSON}.
#' 
#' mongo.bson objects have "mongo.bson" as their class and contain an
#' externally managed pointer to the actual document data. This pointer is
#' stored in the "mongo.bson" attribute of the object.
#' 
#' Objects of class "\link{mongo.bson.iterator}" are used to iterate over a
#' mongo.bson object to enumerate its keys and values.
#' 
#' Objects of class "\link{mongo.bson.buffer}" are used to build BSON
#' documents.
#' 
#' 
#' @name mongo.bson
#' @docType class
#' @seealso \code{\link{mongo.bson.from.list}},\cr
#' \code{\link{mongo.bson.to.list}},\cr \link{mongo.bson.iterator},\cr
#' \link{mongo.bson.buffer},\cr \code{\link{mongo.bson.from.buffer}},\cr
#' \code{\link{mongo.bson.empty}},\cr \code{\link{mongo.find.one}},\cr
#' \code{\link{mongo.bson.destroy}}, \code{link{mongo.shorthand}}.
#' @examples
#' 
#' b <- mongo.bson.from.list(list(name="Fred", age=29, city="Boston"))
#' iter <- mongo.bson.iterator.create(b)  # b is of class "mongo.bson"
#' while (mongo.bson.iterator.next(iter))
#'     print(mongo.bson.iterator.value(iter))
#' 
NULL





#' The mongo.code class
#' 
#' Objects of class "mongo.code" are used to represent javascript code values
#' in BSON documents.
#' 
#' mongo.code objects' value is a string representing the value of the code.
#' 
#' mongo.code objects have "mongo.code" as their class so that\cr
#' \code{\link{mongo.bson.buffer.append}()} may detect them and append the
#' appropriate BSON code-typed value to a buffer.
#' 
#' These mongo.code values may also be present in a list and\cr will be handled
#' properly by \code{\link{mongo.bson.buffer.append.list}()} and\cr
#' \code{\link{mongo.bson.from.list}()}.
#' 
#' 
#' @name mongo.code
#' @docType class
#' @seealso \code{\link{mongo.code.create}},\cr
#' \code{\link{mongo.bson.buffer.append}},\cr
#' \code{\link{mongo.bson.buffer.append.list}},\cr \link{mongo.bson.buffer},\cr
#' \link{mongo.bson}.
#' @examples
#' 
#' buf <- mongo.bson.buffer.create()
#' code <- mongo.code.create("y = x")
#' mongo.bson.buffer.append(buf, "Code", code)
#' lst <- list(c1 = code, One = 1)
#' mongo.bson.buffer.append.list(buf, "listWcode", lst)
#' mongo.bson.buffer.append.code(buf, "Code2", "a = 1")
#' b <- mongo.bson.from.buffer(buf)
#' 
#' # the above will create a mongo.bson object of the following form:
#' # { "Code": (CODE) "y = x", 
#' #   "listWcode" : { "c1"  : (CODE) "y = x",
#' #                   "One" : 1 },
#' #   "Code2" : (CODE) "a = 1" }
#' 
NULL





#' The mongo.code.w.scope class
#' 
#' Objects of class "mongo.code.w.scope" are used to represent javascript code
#' values with scopes in BSON documents.
#' 
#' mongo.code.w.scope objects' value is a string representing the value of the
#' code.
#' 
#' The scope is a \link{mongo.bson} object and is stored in the "scope"
#' attribute of the mongo.code.w.scope object.
#' 
#' mongo.code.w.scope objects have "mongo.code.w.scope" as their class so
#' that\cr \code{\link{mongo.bson.buffer.append}()} may detect them and append
#' the appropriate BSON code-typed value and scope to a buffer.
#' 
#' These mongo.code.w.scope values may also be present in a list and will be
#' handled properly by \code{\link{mongo.bson.buffer.append.list}()} and
#' \code{\link{mongo.bson.from.list}()}.
#' 
#' 
#' @name mongo.code.w.scope
#' @docType class
#' @seealso \code{\link{mongo.code.w.scope.create}},\cr
#' \code{\link{mongo.bson.buffer.append}},\cr
#' \code{\link{mongo.bson.buffer.append.list}},\cr \link{mongo.bson.buffer},\cr
#' \link{mongo.bson}.
#' @examples
#' 
#' buf <- mongo.bson.buffer.create()
#' mongo.bson.buffer.append(buf, "sv", "sx")
#' scope <- mongo.bson.from.buffer(buf)
#' codeWscope <- mongo.code.w.scope.create("y = x", scope)
#' buf <- mongo.bson.buffer.create()
#' mongo.bson.buffer.append(buf, "CodeWscope", codeWscope)
#' lst <- list(c1 = codeWscope, One = 1)
#' mongo.bson.buffer.append.list(buf, "listWcodeWscope", lst)
#' b <- mongo.bson.from.buffer(buf)
#' 
#' # the above will create a mongo.bson object of the following form:
#' # { "CodeWscope" : (CODEWSCOPE) "y = x"
#' #                  (SCOPE) { "sv" : "sx"}, 
#' #   "listWcodeWscope" : { "c1" : (CODEWSCOPE) "y = x"
#' #                                (SCOPE) { "sv" : "sx"} } }
#' 
NULL





#' The mongo.cursor class
#' 
#' Objects of class "mongo.cursor" are returned from \code{\link{mongo.find}()}
#' and used to iterate over the records matching the query.
#' 
#' \code{\link{mongo.cursor.next}(cursor)} is used to step to the first or next
#' record.
#' 
#' \code{\link{mongo.cursor.value}(cursor)} returns a mongo.bson object
#' representing the current record.
#' 
#' \code{\link{mongo.cursor.destroy}(cursor)} releases the resources attached
#' to the cursor.
#' 
#' mongo.cursor objects have "mongo.cursor" as their class and contain an
#' externally managed pointer to the actual cursor data. This pointer is stored
#' in the "mongo.cursor" attribute of the object.
#' 
#' 
#' @name mongo.cursor
#' @docType class
#' @seealso \code{\link{mongo.find}},\cr \code{\link{mongo.cursor.next}},\cr
#' \code{\link{mongo.cursor.value}},\cr \code{\link{mongo.cursor.destroy}}.
#' @examples
#' 
#' mongo <- mongo.create()
#' if (mongo.is.connected(mongo)) {
#'     buf <- mongo.bson.buffer.create()
#'     mongo.bson.buffer.append(buf, "city", "St. Louis")
#'     query <- mongo.bson.from.buffer(buf)
#' 
#'     # Find the first 1000 records in collection people
#'     # of database test where city == "St. Louis"
#'     cursor <- mongo.find(mongo, "test.people", query, limit=1000L)
#'     # Step though the matching records and display them
#'     while (mongo.cursor.next(cursor))
#'         print(mongo.cursor.value(cursor))
#'     mongo.cursor.destroy(cursor)
#' }
#' 
NULL





#' The mongo.gridfile class
#' 
#' Objects of class "mongo.gridfile" are used to access gridfiles on a MongoDB
#' server.  They are created by \code{\link{mongo.gridfs.find}()}.
#' 
#' mongo.gridfile objects have "mongo.gridfile" as their class and contain an
#' externally managed pointer to the actual data used to manage operations on
#' the gridfile. This pointer is stored in the "mongo.gridfile" attribute of
#' the object.  The object also has a "mongo.gridfs" attribute holding a
#' pointer to the mongo.gridfs object used in creation to prevent garbage
#' collection on the mongo.gridfs object while the mongo.gridfile is still
#' active.
#' 
#' 
#' @name mongo.gridfile
#' @docType class
#' @seealso \code{\link{mongo.gridfs}},\cr \code{\link{mongo.gridfs.find}},\cr
#' \code{\link{mongo.gridfile.get.filename}},\cr
#' \code{\link{mongo.gridfile.get.descriptor}},\cr
#' \code{\link{mongo.gridfile.get.filename}},\cr
#' \code{\link{mongo.gridfile.get.length}},\cr
#' \code{\link{mongo.gridfile.get.chunk.size}},\cr
#' \code{\link{mongo.gridfile.get.chunk.count}},\cr
#' \code{\link{mongo.gridfile.get.content.type}},\cr
#' \code{\link{mongo.gridfile.get.upload.date}},\cr
#' \code{\link{mongo.gridfile.get.md5}},\cr
#' \code{\link{mongo.gridfile.get.metadata}},\cr
#' \code{\link{mongo.gridfile.get.chunk}},\cr
#' \code{\link{mongo.gridfile.get.chunks}},\cr
#' \code{\link{mongo.gridfile.read}},\cr \code{\link{mongo.gridfile.seek}},\cr
#' \code{\link{mongo.gridfile.pipe}}.
#' @examples
#' 
#' mongo <- mongo.create()
#' if (mongo.is.connected(mongo)) {
#'     gridfs <- mongo.gridfs.create(mongo, "grid")
#'     mongo.gridfs.store.file(gridfs, "tests/test.R", "test.R")
#' 
#'     gf <- mongo.gridfs.find(gridfs, "test.R")
#'     if( !is.null(gf)){
#'         gf
#'         mongo.gridfile.destroy(gf)
#'     }
#'     mongo.gridfs.destroy(gridfs)
#' }
#' 
NULL





#' The mongo.gridfile.writer class
#' 
#' Objects of class "mongo.gridfile.writer" are used to buffer multiple writes
#' to a single GridFS file.
#' 
#' Use \code{\link{mongo.gridfile.writer.create}} to create an object of this
#' class,\cr \code{\link{mongo.gridfile.writer.write}} to write data to it,
#' and\cr \code{\link{mongo.gridfile.writer.finish}} when done writing.
#' 
#' mongo.gridfile.writer objects have "mongo.gridfile.writer" as their class
#' and contain an externally managed pointer to the actual data used to manage
#' operations on the GridFS. This pointer is stored in the "mongo.gridfile"
#' attribute of the object. The object also has a "mongo.gridfs" attribute
#' holding a pointer to the mongo.gridfs object used in creation to prevent
#' garbage collection on the mongo.gridfs object while the
#' mongo.gridfile.writer is still active.
#' 
#' 
#' @name mongo.gridfile.writer
#' @docType class
#' @seealso \link{mongo.gridfs},\cr
#' \code{\link{mongo.gridfile.writer.create}},\cr
#' \code{\link{mongo.gridfile.writer.write}},\cr
#' \code{\link{mongo.gridfile.writer.finish}}.
#' @examples
#' 
#' mongo <- mongo.create()
#' if (mongo.is.connected(mongo)) {
#'     gridfs <- mongo.gridfs.create(mongo, "grid")
#' 
#'     gfw <- mongo.gridfile.writer.create(gridfs, "test.dat")
#' 
#'     # store 4 bytes
#'     mongo.gridfile.writer.write(gfw, charToRaw("test"))
#' 
#'     # store string & LF plus 0-byte terminator
#'     buf <- writeBin("Test\n", as.raw(1))
#'     mongo.gridfile.writer.write(gfw, buf)
#' 
#'     # store PI as a float
#'     buf <- writeBin(3.1415926, as.raw(1), size=4, endian="little")
#'     mongo.gridfile.writer.write(gfw, buf)
#' 
#'     mongo.gridfile.writer.finish(gfw)
#'     mongo.gridfs.destroy(gridfs)
#' }
#' 
NULL





#' The mongo.gridfs class
#' 
#' Objects of class "mongo.gridfs" are used to store and/or access a "Grid File
#' System" (GridFS) on a MongoDB server.  While primarily intended to store
#' large documents that won't fit on the server as a single BSON object, GridFS
#' may also be used to store large numbers of smaller files.
#' 
#' See \url{http://www.mongodb.org/display/DOCS/GridFS} and\cr
#' \url{http://www.mongodb.org/display/DOCS/When+to+use+GridFS}.
#' 
#' mongo.gridfs objects have "mongo.gridfs" as their class and contain an
#' externally managed pointer to the actual data used to manage operations on
#' the GridFS.\cr This pointer is stored in the "mongo.gridfs" attribute of the
#' object.  The object also has a "mongo" attribute holding a pointer to the
#' mongo connection object used in creation to prevent garbage collection on
#' the mongo object while the mongo.gridfile is still active.
#' 
#' Objects of class "\link{mongo.gridfile}" are used to access gridfiles and
#' read from them.
#' 
#' Objects of class "\link{mongo.gridfile.writer}" are used to write data to
#' the GridFS.
#' 
#' 
#' @name mongo.gridfs
#' @docType class
#' @seealso \code{\link{mongo.gridfs.destroy}},\cr
#' \code{\link{mongo.gridfs.store.file}},\cr
#' \code{\link{mongo.gridfs.remove.file}},\cr
#' \code{\link{mongo.gridfs.store}},\cr
#' \code{\link{mongo.gridfile.writer.create}},\cr
#' \code{\link{mongo.gridfs.find}}.
#' @examples
#' 
#' mongo <- mongo.create()
#' if (mongo.is.connected(mongo)) {
#'     gridfs <- mongo.gridfs.create(mongo, "grid")
#'     # Copy a local file to the server as a GridFS file
#'     mongo.gridfs.store.file(gridfs, "tests/test.R", "test.R")
#' 
#'     # locate the file on the server
#'     gf <- mongo.gridfs.find(gridfs, "test.R")
#'     if( !is.null(gf)){
#'       print(gf)
#'       # and pipe it to an R connection object
#'       test.R <- file("test2.R")
#'       mongo.gridfile.pipe(gf, test.R)
#'   
#'       mongo.gridfile.destroy(gf)
#'     }
#'     mongo.gridfs.destroy(gridfs)
#' }
#' 
NULL





#' The mongo.oid class
#' 
#' Objects of class "mongo.oid" represent MongoDB Object IDs.
#' 
#' See \url{http://www.mongodb.org/display/DOCS/Object+IDs}
#' 
#' mongo.oid objects contain an externally managed pointer to the actual
#' 12-byte object ID data. This pointer is stored in the "mongo.oid" attribute
#' of the object.
#' 
#' mongo.oid objects have "mongo.oid" as their class so that
#' \code{\link{mongo.bson.buffer.append}()} may detect them and append the
#' appropriate BSON OID-typed value to a buffer.
#' 
#' mongo.oid values may also be present in a list and will be handled
#' properly\cr by \code{\link{mongo.bson.buffer.append.list}()} and
#' \code{\link{mongo.bson.from.list}()}.
#' 
#' 
#' @name mongo.oid
#' @docType class
#' @seealso \link{mongo.oid},\cr \code{\link{mongo.oid.from.string}},\cr
#' \code{\link{as.character.mongo.oid}},\cr
#' \code{\link{mongo.oid.to.string}},\cr \code{\link{mongo.oid.time}},\cr
#' \code{\link{mongo.bson.buffer.append}},\cr
#' \code{\link{mongo.bson.buffer.append.oid}},\cr
#' \code{\link{mongo.bson.buffer.append.list}},\cr \link{mongo.bson.buffer},\cr
#' \link{mongo.bson}.
#' @examples
#' 
#' buf <- mongo.bson.buffer.create()
#' oid <- mongo.oid.create()
#' mongo.bson.buffer.append(buf, "_id", oid)
#' b <- mongo.bson.from.buffer(buf)
#' 
NULL





#' The mongo.regex class
#' 
#' Objects of class "mongo.regex" represent regular expressions and are strings
#' with the options value stored in the "options" attribute.
#' 
#' See
#' \url{http://www.mongodb.org/display/DOCS/Advanced+Queries#AdvancedQueries-RegularExpressions}
#' 
#' mongo.regex objects have "mongo.regex" as their class so that\cr
#' \code{\link{mongo.bson.buffer.append}()} may detect them and append the
#' appropriate BSON regex-typed value to a buffer.
#' 
#' These mongo.regex values may also be present in a list and will be handled
#' properly\cr by \code{\link{mongo.bson.buffer.append.list}()} and
#' \code{\link{mongo.bson.from.list}()}.
#' 
#' 
#' @name mongo.regex
#' @docType class
#' @seealso \code{\link{mongo.regex.create}},\cr
#' \code{\link{mongo.bson.buffer.append}},\cr
#' \code{\link{mongo.bson.buffer.append.list}},\cr \link{mongo.bson.buffer},\cr
#' \link{mongo.bson}.
#' @examples
#' 
#' buf <- mongo.bson.buffer.create()
#' regex <- mongo.regex.create("acme.*corp", options="i")
#' mongo.bson.buffer.append.regex(buf, "MatchAcme", regex)
#' b <- mongo.bson.from.buffer(buf)
#' 
NULL





#' The mongo.symbol class
#' 
#' Objects of class "mongo.symbol" are used to represent symbol values in BSON
#' documents.
#' 
#' mongo.symbol objects' value is a string representing the value of the
#' symbol.
#' 
#' mongo.symbol objects have "mongo.symbol" as their class so that\cr
#' \code{\link{mongo.bson.buffer.append}()} may detect them and append the
#' appropriate BSON symbol-typed value to a buffer.
#' 
#' These mongo.symbol values may also be present in a list and will be handled
#' properly\cr by \code{\link{mongo.bson.buffer.append.list}()} and
#' \code{\link{mongo.bson.from.list}()}.
#' 
#' 
#' @name mongo.symbol
#' @docType class
#' @seealso \code{\link{mongo.symbol.create}},\cr
#' \code{\link{mongo.bson.buffer.append}},\cr
#' \code{\link{mongo.bson.buffer.append.list}},\cr \link{mongo.bson.buffer},\cr
#' \link{mongo.bson}.
#' @examples
#' 
#' buf <- mongo.bson.buffer.create()
#' sym <- mongo.symbol.create("Beta")
#' mongo.bson.buffer.append(buf, "B", sym)
#' l <- list(s1 = sym, Two = 2)
#' mongo.bson.buffer.append.list(buf, "listWsym", l)
#' b <- mongo.bson.from.buffer(buf)
#' 
#' # the above will create a mongo.bson object of the following form:
#' # { "B": (SYMBOL) "Beta", 
#' #   "listWsym" : { "s1" : (SYMBOL) "Beta",
#' #                  "Two" : 2 } }
#' 
NULL





#' The mongo.timestamp class
#' 
#' Objects of class "mongo.timestamp" are an extension of the POSIXct class.
#' They have their increment value stored in the "increment" attribute of the
#' object.
#' 
#' See \url{http://www.mongodb.org/display/DOCS/Timestamp+Data+Type}
#' 
#' mongo.timestamp objects have "mongo.timestamp", "POSIXct" & "POSIXt" as
#' their class so that \code{\link{mongo.bson.buffer.append}()} may detect them
#' and append the appropriate BSON code-typed value to a buffer.
#' 
#' These mongo.timestamp values may also be present in a list and will be
#' handled properly by \code{\link{mongo.bson.buffer.append.list}()} and
#' \code{\link{mongo.bson.from.list}()}.
#' 
#' 
#' @name mongo.timestamp
#' @docType class
#' @seealso \code{\link{mongo.timestamp.create}},\cr
#' \code{\link{mongo.bson.buffer.append}},\cr
#' \code{\link{mongo.bson.buffer.append.list}},\cr \link{mongo.bson.buffer},\cr
#' \link{mongo.bson}.
#' @examples
#' 
#' mongo <- mongo.create()
#' if (mongo.is.connected(mongo)) {
#'     buf <- mongo.bson.buffer.create()
#'     # special Null timestamp -- automatically filled in 
#'     # if one of first two fields in a record
#'     ts <- mongo.timestamp.create(0,0)
#'     mongo.bson.buffer.append(buf, "InsertTime", ts)
#'     mongo.bson.buffer.append(buf, "name", "Joe")
#'     b <- mongo.bson.from.buffer(buf)
#'     mongo.insert(mongo, "test.people", b)
#' 
#'     # create using a POSIXlt
#'     ts <- mongo.timestamp.create(strptime("05-12-2012",
#'         "%m-%d-%Y"), increment=1)
#' }
#' 
NULL





#' The mongo.undefined class
#' 
#' Objects of class "mongo.undefined" are used to represent undefined values in
#' BSON documents.
#' 
#' mongo.undefined objects are strings (a character vector) with a single value
#' of "UNDEFINED"
#' 
#' mongo.undefined objects have "mongo.undefined" as their class so that\cr
#' \code{\link{mongo.bson.buffer.append}()} may detect them and append the
#' appropriate BSON undefined value to a buffer.
#' 
#' These mongo.undefined values may also be present in a list and will be
#' handled properly by \code{\link{mongo.bson.buffer.append.list}()} and
#' \code{\link{mongo.bson.from.list}()}.
#' 
#' 
#' @name mongo.undefined
#' @docType class
#' @seealso \code{\link{mongo.undefined.create}},\cr
#' \code{\link{mongo.bson.buffer.append}},\cr
#' \code{\link{mongo.bson.buffer.append.list}},\cr \link{mongo.bson.buffer},\cr
#' \link{mongo.bson}.
#' @examples
#' 
#' buf <- mongo.bson.buffer.create()
#' undef <- mongo.undefined.create()
#' mongo.bson.buffer.append(buf, "Undef", undef)
#' l <- list(u1 = undef, One = 1)
#' mongo.bson.buffer.append.list(buf, "listWundef", l)
#' b <- mongo.bson.from.buffer(buf)
#' 
#' # the above will create a mongo.bson object of the following form:
#' # { "Undef": UNDEFINED, "listWundef" : { "u1" : UNDEFINED, "One" : 1 } }
#' 
NULL





#' The mongo (database connection) class
#' 
#' Objects of class "mongo" are used to connect to a MongoDB server and to
#' perform database operations on that server.
#' 
#' mongo objects have "mongo" as their class and contain an externally managed
#' pointer to the connection data. This pointer is stored in the "mongo"
#' attribute of the object.
#' 
#' Note that the members of the mongo object only reflect\cr the initial
#' parameters of \code{\link{mongo.create}()}. Only the external data actually
#' changes if, for example, mongo.timeout is called after the initial call to
#' \code{mongo.create}.
#' 
#' 
#' @name mongo
#' @docType class
#' @seealso \code{\link{mongo.create}},\cr \code{\link{mongo.is.connected}},\cr
#' \code{\link{mongo.get.databases}},\cr
#' \code{\link{mongo.get.database.collections}},\cr
#' \code{\link{mongo.insert}},\cr \code{\link{mongo.find.one}},\cr
#' \code{\link{mongo.find}},\cr \code{\link{mongo.update}},\cr
#' \code{\link{mongo.remove}},\cr \code{\link{mongo.drop}},\cr
#' \code{\link{mongo.drop.database}}\cr \link{mongo.gridfs}.
#' @examples
#' 
#' mongo <- mongo.create()
#' if (mongo.is.connected(mongo)) {
#'     buf <- mongo.bson.buffer.create()
#'     mongo.bson.buffer.append(buf, "name", "Joe")
#'     mongo.bson.buffer.append(buf, "age", 22L)
#'     b <- mongo.bson.from.buffer(buf)
#'     mongo.insert(mongo, "test.people", b)
#' }
#' 
NULL






#' zips Dataset
#' 
#' A dataset with US zip data provided by mongodb for education.
#' 
#' 
#' @name zips
#' @aliases zips Dataset zips
#' @docType data
#' @format The _id field holds the zip code as a string.
#' 
#' The city field holds the city.
#' 
#' The state field holds the two letter state abbreviation.
#' 
#' The pop field holds the population.
#' 
#' The loc field holds the location as a latitude longitude pair.
#' @source \url{http://media.mongodb.org/zips.json}
#' @examples
#' 
#' \dontrun{
#' # code to create the dataset
#' library(RJSONIO)
#' 
#' json_file <- "http://media.mongodb.org/zips.json"
#' 
#' rL <- readLines(json_file)
#' 
#' zips <- do.call(rbind,lapply(rL,fromJSON))
#' 
#' save(zips, file="data/zips.rda", compress="xz")
#' }
#' 
NULL


