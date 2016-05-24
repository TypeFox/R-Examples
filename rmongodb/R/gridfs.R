
#' Create a mongo.gridfs object
#' 
#' Create a \link{mongo.gridfs} object used to access and store "grid files" on
#' the MongoDB server.
#' 
#' 
#' @param mongo A (\link{mongo}) connection object.
#' @param db (string) The name of the database in which to access and/or store
#' the gridfs-related collections.
#' @param prefix (string) The prefix to use constructing the gridfs-related
#' collection names.  There are two collections used for this purpose:\cr
#' \"\code{db}.\code{prefix}.files\" and \"\code{db}.\code{prefix}.chunks\".
#' @return (\link{mongo.gridfs}) An object to be used for subsequent operations
#' on the grid file store.
#' @seealso \link{mongo.gridfs},\cr \code{\link{mongo.gridfs.destroy}},\cr
#' \code{\link{mongo.gridfs.store.file}},\cr
#' \code{\link{mongo.gridfs.remove.file}},\cr
#' \code{\link{mongo.gridfs.store}},\cr
#' \code{\link{mongo.gridfile.writer.create}},\cr
#' \code{\link{mongo.gridfs.find}}, \code{link{mongo.shorthand}}.
#' @examples
#' 
#' mongo <- mongo.create()
#' if (mongo.is.connected(mongo)) {
#'     gridfs <- mongo.gridfs.create(mongo, "grid")
#'     # Copy a local file to the server as a gridfs file
#'     mongo.gridfs.store.file(gridfs, "tests/test.R", "test.R")
#'     mongo.gridfs.destroy(gridfs)
#' }
#' 
#' @export mongo.gridfs.create
mongo.gridfs.create <- function(mongo, db, prefix="fs")
    .Call(".mongo.gridfs.create", mongo, db, prefix)



#' Destroy a mongo.gridfs object
#' 
#' Releases the resources associated with a \link{mongo.gridfs} object.
#' 
#' It is not absolutely necessary to call this function since R's garbage
#' collection will eventually get around to doing it for you.
#' 
#' 
#' @param gridfs A (\link{mongo.gridfs}) object.
#' @return NULL
#' @seealso \link{mongo.gridfs},\cr \code{\link{mongo.gridfs.create}},\cr
#' \code{\link{mongo.gridfs.store.file}}.
#' @examples
#' 
#' mongo <- mongo.create()
#' if (mongo.is.connected(mongo)) {
#'     gridfs <- mongo.gridfs.create(mongo, "grid")
#'     # Copy a local file to the server as a gridfs file
#'     mongo.gridfs.store.file(gridfs, "tests/test.R", "test.R")
#'     mongo.gridfs.destroy(gridfs)
#' }
#' 
#' @export mongo.gridfs.destroy
mongo.gridfs.destroy <- function(gridfs)
    .Call(".mongo.gridfs.destroy", gridfs)



#' Store a file into a GridFS on a MongoDB server
#' 
#' Store a file into a GridFS on a MongoDB server.  This function stores the
#' entire given file on the server, breaking it up into 256K chunks as
#' necessary.
#' 
#' 
#' @param gridfs A (\link{mongo.gridfs}) object.
#' @param filename (string) The path/filename of the file to copy to the
#' server.
#' @param remotename (string) The name the file will be known as within the
#' GridFS.\cr If remotename=="" (the default), the remote file will be known by
#' the given \code{filename}.
#' @param contenttype (string) Optional MIME content type.
#' @return TRUE, if successful; FALSE, if an error occured during the
#' operation.
#' @seealso \link{mongo.gridfs},\cr \code{\link{mongo.gridfs.create}},\cr
#' \code{\link{mongo.gridfs.remove.file}}.
#' @examples
#' 
#' mongo <- mongo.create()
#' if (mongo.is.connected(mongo)) {
#'     gridfs <- mongo.gridfs.create(mongo, "grid")
#'     # Copy a local file to the server as a gridfs file
#'     mongo.gridfs.store.file(gridfs, "tests/test.R", "test.R")
#'     mongo.gridfs.destroy(gridfs)
#' }
#' 
#' @export mongo.gridfs.store.file
mongo.gridfs.store.file <- function(gridfs, filename, remotename="", contenttype="")
    .Call(".mongo.gridfs.store.file", gridfs, filename, remotename, contenttype)



#' Remove a file from a GridFS on a MongoDB server
#' 
#' Remove a file from a GridFS on a MongoDB server.
#' 
#' 
#' @param gridfs A (\link{mongo.gridfs}) object.
#' @param remotename (string) The name of the file to be removed (as known
#' within the GridFS).
#' @return NULL
#' @seealso \link{mongo.gridfs},\cr \code{\link{mongo.gridfs.store.file}}\cr
#' \code{\link{mongo.gridfs.store}}.
#' @examples
#' 
#' mongo <- mongo.create()
#' if (mongo.is.connected(mongo)) {
#'     gridfs <- mongo.gridfs.create(mongo, "grid")
#'     \dontrun{mongo.gridfs.remove.file(gridfs, "test.R")}
#'     mongo.gridfs.destroy(gridfs)
#' }
#' 
#' @export mongo.gridfs.remove.file
mongo.gridfs.remove.file <- function(gridfs, remotename)
    .Call(".mongo.gridfs.remove.file", gridfs, remotename)



#' Store raw data as a file in a GridFS
#' 
#' Store raw data as a file to a GridFS on a MongoDB server. This function
#' stores the entire piece of data file on the server, breaking it up into 256K
#' chunks as necessary.
#' 
#' This function only handles the RAW type. Use \code{writeBin()} as necessary
#' to pack your data appropriately for storage.  See the examples and R's
#' documentation on \code{writeBin()}.
#' 
#' Use \link{mongo.gridfile.writer} when you need to buffer many writes to a
#' GridFS file.
#' 
#' 
#' @param gridfs A (\link{mongo.gridfs}) object.
#' @param raw (raw) The data to store on the server.
#' @param remotename (string) The name the file will be known as within the
#' GridFS.
#' @param contenttype (string) Optional MIME content type.
#' @return TRUE, if successful; FALSE, if an error occured during the
#' operation.
#' @seealso \link{mongo.gridfs},\cr \code{\link{mongo.gridfs.create}},\cr
#' \code{\link{mongo.gridfs.remove.file}}.
#' @examples
#' 
#' mongo <- mongo.create()
#' if (mongo.is.connected(mongo)) {
#'     gridfs <- mongo.gridfs.create(mongo, "grid")
#'     # store 4 bytes
#'     mongo.gridfs.store(gridfs, charToRaw("test"), "test4.dat")
#' 
#'     # store string & LF plus 0-byte terminator
#'     buf <- writeBin("Test\n", as.raw(1))
#'     mongo.gridfs.store(gridfs, buf, "test6.dat")
#' 
#'     # store PI as a float
#'     buf <- writeBin(3.1415926, as.raw(1), size=4, endian="little")
#'     mongo.gridfs.store(gridfs, buf, "PI.dat")
#' 
#'     mongo.gridfs.destroy(gridfs)
#' }
#' 
#' @export mongo.gridfs.store
mongo.gridfs.store <- function(gridfs, raw, remotename, contenttype="")
    .Call(".mongo.gridfs.store", gridfs, raw, remotename, contenttype)



#' Create a mongo.gridfile.writer object
#' 
#' Create a \link{mongo.gridfile.writer} object used to buffer many writes to a
#' single GridFS file. Once the mongo.gridfile.writer is created, use
#' \code{\link{mongo.gridfile.writer.write}()} to write data to the buffered
#' GridFS file and \code{\link{mongo.gridfile.writer.finish}()} when done.
#' 
#' 
#' @param gridfs A (\link{mongo.gridfs}) object.
#' @param remotename (string) The name the file will be known as within the
#' GridFS.
#' @param contenttype (string) Optional MIME content type.
#' @return (\link{mongo.gridfile.writer}) The object to be used for writing to
#' the GridFS file.
#' @seealso \link{mongo.gridfs},\cr \code{\link{mongo.gridfs.create}},\cr
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
#' @export mongo.gridfile.writer.create
mongo.gridfile.writer.create <- function(gridfs, remotename, contenttype="")
    .Call(".mongo.gridfile.writer.create", gridfs, remotename, contenttype)



#' Write raw data to a buffered GridFS file
#' 
#' Write raw data to a buffered GridFS file. The data is buffered and sent to
#' the server in 256k chunks as it accumulates.
#' 
#' This function only handles the RAW type. Use \code{writeBin()} as necessary
#' to pack your data appropriately for storage.  See the examples and R's
#' documentation on \code{writeBin()}.
#' 
#' Use \code{\link{mongo.gridfs.store}()} when you only need to write one data
#' packet as a complete GridFS file.
#' 
#' 
#' @param gfw A (\link{mongo.gridfile.writer}) object.
#' @param raw (raw) The data to write to the GridFS file.
#' @return NULL
#' @seealso \link{mongo.gridfs},\cr
#' \code{\link{mongo.gridfile.writer.create}},\cr
#' \link{mongo.gridfile.writer},\cr \code{\link{mongo.gridfile.writer.finish}}.
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
#' 
#'     mongo.gridfs.destroy(gridfs)
#' }
#' 
#' @export mongo.gridfile.writer.write
mongo.gridfile.writer.write <- function(gfw, raw)
    .Call(".mongo.gridfile.writer.write", gfw, raw)



#' Finish writing to a buffered GridFS file
#' 
#' Finish writing to a buffered GridFS file. This function flushes any partial
#' buffer and finalizes the operation.
#' 
#' 
#' @param gfw A (\link{mongo.gridfile.writer}) object.
#' @return TRUE, if successfil; false, if an error occurred.
#' @seealso \link{mongo.gridfs},\cr
#' \code{\link{mongo.gridfile.writer.create}},\cr
#' \link{mongo.gridfile.writer},\cr \code{\link{mongo.gridfile.writer.write}}.
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
#' @export mongo.gridfile.writer.finish
mongo.gridfile.writer.finish <- function(gfw)
    .Call(".mongo.gridfile.writer.finish", gfw)



#' Find a GridFS file
#' 
#' Find a GridFS file and return a \link{mongo.gridfile} object used for
#' further operations on it
#' 
#' 
#' @param gridfs A (\link{mongo.gridfs}) object.
#' @param query (string) The name of the GridFS file to locate.
#' 
#' This parameter may also be a \link{mongo.bson} query object and is used to
#' search the GridFS "files" collection documents for matches.  Alternately,
#' \code{query} may be a list which will be converted to a mongo.bson object by
#' \code{\link{mongo.bson.from.list}()}.
#' @return NULL, if not found; otherwise, a \link{mongo.gridfile} object
#' corresponding to the found GridFS file.
#' @seealso \link{mongo.gridfile},\cr
#' \code{\link{mongo.gridfile.get.filename}},\cr \link{mongo.gridfs}.
#' @examples
#' 
#' mongo <- mongo.create()
#' if (mongo.is.connected(mongo)) {
#'     gridfs <- mongo.gridfs.create(mongo, "grid")
#' 
#'     gf <- mongo.gridfs.find(gridfs, "test.dat")
#'     print(mongo.gridfile.get.length(gf))
#' 
#'     # find a GridFS file uploaded midnight July 4, 2008
#'     buf <- mongo.bson.buffer.create()
#'     mongo.bson.buffer.append(buf, "uploadDate",
#'         strptime("07-04-2008", "%m-%d-%Y"))
#'     query <- mongo.bson.from.buffer(buf)
#'     gf <- mongo.gridfs.find(gridfs, query)
#' 
#'     if (!is.null(gf))
#'         print(mongo.gridfile.get.filename(gf))
#' 
#'     mongo.gridfs.destroy(gridfs)
#' }
#' 
#' @export mongo.gridfs.find
mongo.gridfs.find <- function(gridfs, query) {
    if (typeof(query) == "list")
        query <- mongo.bson.from.list(query)
    .Call(".mongo.gridfs.find", gridfs, query)
}



#' Destroy a mongo.gridfile object
#' 
#' Releases the resources associated with a \link{mongo.gridfile} object.\cr
#' These are created by \code{\link{mongo.gridfs.find}()}.
#' 
#' It is not absolutely necessary to call this function since R's garbage
#' collection will eventually get around to doing it for you.
#' 
#' 
#' @param gridfile A (\link{mongo.gridfile}) object.
#' @return NULL
#' @seealso \code{\link{mongo.gridfs.find}},\cr \link{mongo.gridfile},\cr
#' \link{mongo.gridfs}.
#' @examples
#' 
#' mongo <- mongo.create()
#' if (mongo.is.connected(mongo)) {
#'     gridfs <- mongo.gridfs.create(mongo, "grid")
#'     mongo.gridfs.store.file(gridfs, "tests/test.R", "test.R")
#'     gf <- mongo.gridfs.find(gridfs, "test.R")
#'     if( !is.null(gf) ){
#'       print(mongo.gridfile.get.upload.date(gf))
#'       mongo.gridfile.destroy(gf)
#'     }
#'     mongo.gridfs.destroy(gridfs)
#' }
#' 
#' @export mongo.gridfile.destroy
mongo.gridfile.destroy <- function(gridfile)
    .Call(".mongo.gridfile.destroy", gridfile)



#' Get the descriptor of a mongo.gridfile
#' 
#' Get the descriptor of a \link{mongo.gridfile}.  This descriptor is the
#' document describing the given gridfile as stored on the MongoDB server in
#' the 'files' collection of the GridFS .
#' 
#' See \url{http://www.mongodb.org/display/DOCS/GridFS+Specification}.
#' 
#' 
#' @param gridfile A (\link{mongo.gridfile}) object.
#' @return (\link{mongo.bson}) The descriptor of \code{gridfile}.
#' @seealso \code{\link{mongo.gridfs}},\cr \code{\link{mongo.gridfs.find}},\cr
#' \link{mongo.gridfile},\cr \code{\link{mongo.gridfile.get.filename}},\cr
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
#'       print(mongo.gridfile.get.descriptor(gf))
#' 
#'       mongo.gridfile.destroy(gf)
#'     }
#'     mongo.gridfs.destroy(gridfs)
#' }
#' 
#' @export mongo.gridfile.get.descriptor
mongo.gridfile.get.descriptor <- function(gridfile)
    .Call(".mongo.gridfile.get.descriptor", gridfile)



#' Get the filename of a mongo.gridfile
#' 
#' Get the filename of a \link{mongo.gridfile}. This is the 'remote name' that
#' is used identify the file on the server.
#' 
#' 
#' @param gridfile A (\link{mongo.gridfile}) object.
#' @return (string) The filename (remote name) of \code{gridfile}
#' @seealso \code{\link{mongo.gridfs}},\cr \code{\link{mongo.gridfs.find}},\cr
#' \link{mongo.gridfile},\cr \code{\link{mongo.gridfile.get.descriptor}},\cr
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
#' 
#'     # find a GridFS file uploaded midnight July 4, 2008
#'     buf <- mongo.bson.buffer.create()
#'     mongo.bson.buffer.append(buf, "uploadDate",
#'         strptime("07-04-2008", "%m-%d-%Y"))
#'     query <- mongo.bson.from.buffer(buf)
#' 
#'     gf <- mongo.gridfs.find(gridfs, query)
#'     if (!is.null(gf)) {
#'         print(mongo.gridfile.get.filename(gf))
#' 
#'         mongo.gridfile.destroy(gf)
#'     }
#'     mongo.gridfs.destroy(gridfs)
#' }
#' 
#' @export mongo.gridfile.get.filename
mongo.gridfile.get.filename <- function(gridfile)
    .Call(".mongo.gridfile.get.filename", gridfile)



#' Get the length of a mongo.gridfile
#' 
#' Get the length of a \link{mongo.gridfile}.
#' 
#' 
#' @param gridfile A (\link{mongo.gridfile}) object.
#' @return (double) The length of \code{gridfile}.
#' @seealso \code{\link{mongo.gridfs}},\cr \code{\link{mongo.gridfs.find}},\cr
#' \link{mongo.gridfile},\cr \code{\link{mongo.gridfile.get.descriptor}},\cr
#' \code{\link{mongo.gridfile.get.filename}},\cr
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
#'     if( !is.null(gf) ){
#'       print(mongo.gridfile.get.length(gf))
#' 
#'       mongo.gridfile.destroy(gf)
#'     }
#'     mongo.gridfs.destroy(gridfs)
#' }
#' 
#' @export mongo.gridfile.get.length
mongo.gridfile.get.length <- function(gridfile)
    .Call(".mongo.gridfile.get.length", gridfile)



#' Get the chunk.size of a mongo.gridfile
#' 
#' Get the chunk size of a \link{mongo.gridfile}.  This is the size of the
#' chunks into which file is broken up on the server.
#' 
#' 
#' @param gridfile A (\link{mongo.gridfile}) object.
#' @return (integer) The chunk size of \code{gridfile}.
#' @seealso \code{\link{mongo.gridfs}},\cr \code{\link{mongo.gridfs.find}},\cr
#' \link{mongo.gridfile},\cr \code{\link{mongo.gridfile.get.descriptor}},\cr
#' \code{\link{mongo.gridfile.get.filename}},\cr
#' \code{\link{mongo.gridfile.get.length}},\cr
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
#' 
#'     gf <- mongo.gridfs.find(gridfs, "test.R")
#'     if( !is.null(gf)){
#'       print(mongo.gridfile.get.chunk.size(gf))
#' 
#'       mongo.gridfile.destroy(gf)
#'     }
#'     mongo.gridfs.destroy(gridfs)
#' }
#' 
#' @export mongo.gridfile.get.chunk.size
mongo.gridfile.get.chunk.size <- function(gridfile)
    .Call(".mongo.gridfile.get.chunk.size", gridfile)



#' Get the chunk count of a mongo.gridfile
#' 
#' Get the chunk count of a \link{mongo.gridfile}. This is the number of chunks
#' into which the gridfile is broken up on the server.
#' 
#' 
#' @param gridfile A (\link{mongo.gridfile}) object.
#' @return (integer) The chunk count of \code{gridfile}
#' @seealso \code{\link{mongo.gridfs}},\cr \code{\link{mongo.gridfs.find}},\cr
#' \link{mongo.gridfile},\cr \code{\link{mongo.gridfile.get.descriptor}},\cr
#' \code{\link{mongo.gridfile.get.filename}},\cr
#' \code{\link{mongo.gridfile.get.length}},\cr
#' \code{\link{mongo.gridfile.get.chunk.size}},\cr
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
#'       print(mongo.gridfile.get.chunk.count(gf))
#' 
#'       mongo.gridfile.destroy(gf)
#'     }
#'     mongo.gridfs.destroy(gridfs)
#' }
#' 
#' @export mongo.gridfile.get.chunk.count
mongo.gridfile.get.chunk.count <- function(gridfile)
    .Call(".mongo.gridfile.get.chunk.count", gridfile)



#' Get the content type of a mongo.gridfile
#' 
#' Get the MIME content type of a \link{mongo.gridfile}.
#' 
#' 
#' @param gridfile A (\link{mongo.gridfile}) object.
#' @return (string) The content.type (remote name) of \code{gridfile}. This may
#' be an empty string if no content type is associated with the gridfile.
#' @seealso \code{\link{mongo.gridfs}},\cr \code{\link{mongo.gridfs.find}},\cr
#' \link{mongo.gridfile},\cr \code{\link{mongo.gridfile.get.descriptor}},\cr
#' \code{\link{mongo.gridfile.get.filename}},\cr
#' \code{\link{mongo.gridfile.get.length}},\cr
#' \code{\link{mongo.gridfile.get.chunk.size}},\cr
#' \code{\link{mongo.gridfile.get.chunk.count}},\cr
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
#'       print(mongo.gridfile.get.content.type(gf))
#' 
#'       mongo.gridfile.destroy(gf)
#'     }
#'     mongo.gridfs.destroy(gridfs)
#' }
#' 
#' @export mongo.gridfile.get.content.type
mongo.gridfile.get.content.type <- function(gridfile)
    .Call(".mongo.gridfile.get.content.type", gridfile)



#' Get the upload date of a mongo.gridfile
#' 
#' Get the upload date of a \link{mongo.gridfile}.
#' 
#' 
#' @param gridfile A (\link{mongo.gridfile}) object.
#' @return (POSIXct) The upload date/time of \code{gridfile}.
#' @seealso \code{\link{mongo.gridfs}},\cr \code{\link{mongo.gridfs.find}},\cr
#' \link{mongo.gridfile},\cr \code{\link{mongo.gridfile.get.descriptor}},\cr
#' \code{\link{mongo.gridfile.get.filename}},\cr
#' \code{\link{mongo.gridfile.get.length}},\cr
#' \code{\link{mongo.gridfile.get.chunk.size}},\cr
#' \code{\link{mongo.gridfile.get.chunk.count}},\cr
#' \code{\link{mongo.gridfile.get.content.type}},\cr
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
#'     if( !is.null(gf) ){
#'       print(mongo.gridfile.get.upload.date(gf))
#' 
#'       mongo.gridfile.destroy(gf)
#'     }
#'     mongo.gridfs.destroy(gridfs)
#' }
#' 
#' @export mongo.gridfile.get.upload.date
mongo.gridfile.get.upload.date <- function(gridfile)
    .Call(".mongo.gridfile.get.upload.date", gridfile)



#' Get the MD5 hash of a mongo.gridfile
#' 
#' Get the MD5 hash of a \link{mongo.gridfile}.
#' 
#' 
#' @param gridfile A (\link{mongo.gridfile}) object.
#' @return (string) The MD5 hash (32 hex digits) of \code{gridfile}.
#' @seealso \code{\link{mongo.gridfs}},\cr \code{\link{mongo.gridfs.find}},\cr
#' \link{mongo.gridfile},\cr \code{\link{mongo.gridfile.get.descriptor}},\cr
#' \code{\link{mongo.gridfile.get.filename}},\cr
#' \code{\link{mongo.gridfile.get.length}},\cr
#' \code{\link{mongo.gridfile.get.chunk.size}},\cr
#' \code{\link{mongo.gridfile.get.chunk.count}},\cr
#' \code{\link{mongo.gridfile.get.content.type}},\cr
#' \code{\link{mongo.gridfile.get.upload.date}},\cr
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
#'     if( !is.null(gf) ){
#'       print(mongo.gridfile.get.md5(gf))
#' 
#'       mongo.gridfile.destroy(gf)
#'     }
#'     mongo.gridfs.destroy(gridfs)
#' }
#' 
#' @export mongo.gridfile.get.md5
mongo.gridfile.get.md5 <- function(gridfile)
    .Call(".mongo.gridfile.get.md5", gridfile)



#' Get the metadata of a mongo.gridfile
#' 
#' Get the metadata of a \link{mongo.gridfile}.  Some applications may store
#' metadata pertaining to a GridFS file in the "metadata" field of the
#' descriptor.  (See \code{\link{mongo.gridfile.get.descriptor}()}.  This
#' function retrieves that field as a \link{mongo.bson} object.
#' 
#' 
#' @param gridfile A (\link{mongo.gridfile}) object.
#' @return (\link{mongo.bson}) The metadata of \code{gridfile} if present;
#' otherwise, NULL.
#' @seealso \code{\link{mongo.gridfs}},\cr \code{\link{mongo.gridfs.find}},\cr
#' \link{mongo.gridfile},\cr \code{\link{mongo.gridfile.get.descriptor}},\cr
#' \code{\link{mongo.gridfile.get.filename}},\cr
#' \code{\link{mongo.gridfile.get.length}},\cr
#' \code{\link{mongo.gridfile.get.chunk.size}},\cr
#' \code{\link{mongo.gridfile.get.chunk.count}},\cr
#' \code{\link{mongo.gridfile.get.content.type}},\cr
#' \code{\link{mongo.gridfile.get.upload.date}},\cr
#' \code{\link{mongo.gridfile.get.md5}},\cr
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
#'     if( !is.null(gf) ){
#'       print(mongo.gridfile.get.metadata(gf))
#'   
#'       mongo.gridfile.destroy(gf)
#'     }
#'     mongo.gridfs.destroy(gridfs)
#' }
#' 
#' @export mongo.gridfile.get.metadata
mongo.gridfile.get.metadata <- function(gridfile)
    .Call(".mongo.gridfile.get.metadata", gridfile)



#' Get a chunk of a mongo.gridfile
#' 
#' Get a chunk of a \link{mongo.gridfile}.
#' 
#' 
#' @param gridfile A (\link{mongo.gridfile}) object.
#' @param i (integer) The index of the chunk to fetch.  This should be in the
#' range 0 to \code{\link{mongo.gridfile.get.chunk.count}(gridfile) - 1}.
#' @return (\link{mongo.bson}) the \code{i}th chunk of \code{gridfile} if
#' successful; otherwise, NULL.
#' 
#' The value returned is the \code{i}th document in the 'chunks' collection of
#' the GridFS.  The 'data' field of this document contains the actual data
#' belonging to the chunk.\cr See
#' \url{http://www.mongodb.org/display/DOCS/GridFS+Specification}.
#' @seealso \code{\link{mongo.gridfs}},\cr \code{\link{mongo.gridfs.find}},\cr
#' \link{mongo.gridfile},\cr \code{\link{mongo.gridfile.get.descriptor}},\cr
#' \code{\link{mongo.gridfile.get.filename}},\cr
#' \code{\link{mongo.gridfile.get.length}},\cr
#' \code{\link{mongo.gridfile.get.chunk.size}},\cr
#' \code{\link{mongo.gridfile.get.chunk.count}},\cr
#' \code{\link{mongo.gridfile.get.content.type}},\cr
#' \code{\link{mongo.gridfile.get.upload.date}},\cr
#' \code{\link{mongo.gridfile.get.md5}},\cr
#' \code{\link{mongo.gridfile.get.metadata}},\cr
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
#'       chunk <- mongo.gridfile.get.chunk(gf, 0)
#'       iter <- mongo.bson.find(chunk, "data")
#'     
#' 
#'       f <- file("testChunk0.R", "wb")
#'       # write the binary (raw) data to a file
#'       writeBin(mongo.bson.iterator.value(iter), f)
#'       close(f)
#' 
#'       mongo.gridfile.destroy(gf)
#'     }
#'     mongo.gridfs.destroy(gridfs)
#' }
#' 
#' @export mongo.gridfile.get.chunk
mongo.gridfile.get.chunk <- function(gridfile, i)
    .Call(".mongo.gridfile.get.chunk", gridfile, i)



#' Get a cursor for a range of chunks in a mongo.gridfile
#' 
#' Get a \link{mongo.cursor} for a range of chunks in a \link{mongo.gridfile}.
#' 
#' 
#' @param gridfile A (\link{mongo.gridfile}) object.
#' @param start (integer) The index of the first chunk to fetch.  This should
#' be in the range 0 to \code{\link{mongo.gridfile.get.chunk.count}(gridfile) -
#' 1}.
#' @param count (integer) The number of chunks to fetch.
#' @return (\link{mongo.cursor}) A cursor to be used to step through the
#' requested chunks.
#' 
#' The values returned by \code{\link{mongo.cursor.value}()} will be
#' consecutive documents in the 'chunks' collection of the GridFS.  The 'data'
#' field of these documents contains the actual data belonging to the chunk.
#' See \url{http://www.mongodb.org/display/DOCS/GridFS+Specification}.
#' @seealso \code{\link{mongo.gridfs}},\cr \code{\link{mongo.gridfs.find}},\cr
#' \link{mongo.gridfile},\cr \code{\link{mongo.gridfile.get.descriptor}},\cr
#' \code{\link{mongo.gridfile.get.filename}},\cr
#' \code{\link{mongo.gridfile.get.length}},\cr
#' \code{\link{mongo.gridfile.get.chunk.size}},\cr
#' \code{\link{mongo.gridfile.get.chunk.count}},\cr
#' \code{\link{mongo.gridfile.get.content.type}},\cr
#' \code{\link{mongo.gridfile.get.upload.date}},\cr
#' \code{\link{mongo.gridfile.get.md5}},\cr
#' \code{\link{mongo.gridfile.get.metadata}},\cr
#' \code{\link{mongo.gridfile.get.chunk}},\cr
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
#'       cursor <- mongo.gridfile.get.chunks(gf, 1, 2)
#'   
#'       f <- file("rmongodb.pdf.chunks12", "wb")
#'       while (mongo.cursor.next(cursor)) {
#'           chunk <- mongo.cursor.value(cursor)
#'           iter <- mongo.bson.find(chunk, "data")
#'   
#'           # write the binary (raw) data to the file
#'           writeBin(mongo.bson.iterator.value(iter), f)
#'       }
#'       close(f)
#'       mongo.gridfile.destroy(gf)
#'     }
#'     mongo.gridfs.destroy(gridfs)
#' }
#' 
#' @export mongo.gridfile.get.chunks
mongo.gridfile.get.chunks <- function(gridfile, start, count)
    .Call(".mongo.gridfile.get.chunks", gridfile, start, count)



#' Read raw data from a mongo.gridfile
#' 
#' Read raw data from a \link{mongo.gridfile}.  The data read may span multiple
#' chunks.
#' 
#' A mongo.gridfile file maintains a current read position which is advanced by
#' the size of each read.  This position is initially at offset 0.
#' 
#' Since this function returns raw data, you may want to use R's
#' \code{readBin()} to unpack it.
#' 
#' 
#' @param gridfile A (\link{mongo.gridfile}) object.
#' @param size (as.double) The number of bytes to read.
#' @return (raw) The data read from emphgridfile. The length of this vector may
#' be less than the requested size if there was not enough data remaining to be
#' read. This length could also be 0 if an error occured during the operation.
#' Check \code{\link{mongo.get.err}()} of the associated mongo connection
#' object in this case.
#' @seealso \code{\link{mongo.gridfs}},\cr \code{\link{mongo.gridfs.find}},\cr
#' \link{mongo.gridfile},\cr \code{\link{mongo.gridfile.get.descriptor}},\cr
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
#' \code{\link{mongo.gridfile.seek}},\cr \code{\link{mongo.gridfile.pipe}}.
#' @examples
#' 
#' mongo <- mongo.create()
#' if (mongo.is.connected(mongo)) {
#'     gridfs <- mongo.gridfs.create(mongo, "grid")
#'     mongo.gridfs.store.file(gridfs, "tests/test.R", "test.R")
#' 
#'     gf <- mongo.gridfs.find(gridfs, "test.R")
#'     if( !is.null(gf)){
#'       mongo.gridfile.seek(gf, 256*256*5)
#'       data <- mongo.gridfile.read(gf, 16384)
#' 
#'       mongo.gridfile.destroy(gf)
#'     }
#'     mongo.gridfs.destroy(gridfs)
#' }
#' 
#' @export mongo.gridfile.read
mongo.gridfile.read <- function(gridfile, size)
    .Call(".mongo.gridfile.read", gridfile, size)



#' Seek to a position in a mongo.gridfile
#' 
#' Seek to a position in a \link{mongo.gridfile}.\cr This sets the position at
#' which the next \code{\link{mongo.gridfile.read}()} will start.
#' 
#' 
#' @param gridfile A (\link{mongo.gridfile}) object.
#' @param offset (as.double) The position to which to seek.
#' @return (double) The position set.  This may be at the length of the GridFS
#' file if \code{offset} was greater than that.
#' @seealso \code{\link{mongo.gridfs}},\cr \code{\link{mongo.gridfs.find}},\cr
#' \link{mongo.gridfile},\cr \code{\link{mongo.gridfile.get.descriptor}},\cr
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
#' \code{\link{mongo.gridfile.read}},\cr \code{\link{mongo.gridfile.pipe}}.
#' @examples
#' 
#' mongo <- mongo.create()
#' if (mongo.is.connected(mongo)) {
#'     gridfs <- mongo.gridfs.create(mongo, "grid")
#'     mongo.gridfs.store.file(gridfs, "tests/test.R", "test.R")
#' 
#'     gf <- mongo.gridfs.find(gridfs, "test.R")
#'     if( !is.null(gf)){
#'       mongo.gridfile.seek(gf, 256*256*5)
#'       data <- mongo.gridfile.read(gf, 16384)
#'   
#'       mongo.gridfile.destroy(gf)
#'     }
#'     mongo.gridfs.destroy(gridfs)
#' }
#' 
#' @export mongo.gridfile.seek
mongo.gridfile.seek <- function(gridfile, offset)
    .Call(".mongo.gridfile.get.chunk", gridfile, offset)



#' Pipe a mongo.gridfile to an R connection
#' 
#' Pipe a mongo.gridfile to an R connection.  This outputs the entire GridFS
#' file to a connection. If the connection is open, it must be in binary output
#' mode; otherwise, the connection is opened in binary output mode and closed
#' afterwards.
#' 
#' 
#' @param gridfile A (\link{mongo.gridfile}) object.
#' @param con (connection) An R connection object.
#' @return NULL
#' @seealso \code{\link{mongo.gridfs}},\cr \code{\link{mongo.gridfs.find}},\cr
#' \link{mongo.gridfile},\cr \code{\link{mongo.gridfile.get.descriptor}},\cr
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
#' \code{\link{mongo.gridfile.read}},\cr \code{\link{mongo.gridfile.seek}}.
#' @examples
#' 
#' mongo <- mongo.create()
#' if (mongo.is.connected(mongo)) {
#'     gridfs <- mongo.gridfs.create(mongo, "grid")
#'     mongo.gridfs.store.file(gridfs, "tests/test.R", "test.R")
#' 
#'     gf <- mongo.gridfs.find(gridfs, "test.R")
#'     if (!is.null(gf)) {
#'         f <- file("mongodb_copy.pdf")
#'         mongo.gridfile.pipe(gf, f)
#'         
#'         mongo.gridfile.destroy(gf)
#'     }
#' 
#'     mongo.gridfs.destroy(gridfs)
#' }
#' 
#' @export mongo.gridfile.pipe
mongo.gridfile.pipe <- function(gridfile, con) {
    wasOpen <- isOpen(con)
    if (!wasOpen)
        open(con, "wb")
    count <- mongo.gridfile.get.chunk.count(gridfile)
    cursor <- mongo.gridfile.get.chunks(gridfile, 0, count)
    while (mongo.cursor.next(cursor)) {
        b <- mongo.cursor.value(cursor);
        data <- mongo.bson.value(b, "data")
        writeBin(data, con)
        mongo.bson.destroy(b)
    }
    if (!wasOpen)
        close(con)
    mongo.cursor.destroy(cursor)
    NULL
}

