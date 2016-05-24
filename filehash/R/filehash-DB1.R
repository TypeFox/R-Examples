######################################################################
## Copyright (C) 2006--2008, Roger D. Peng <rpeng@jhsph.edu>
##     
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301, USA
#####################################################################

######################################################################
## Class 'filehashDB1'

## Database entries
##
## File format: [key]        [nbytes data] [data]
##              serialized   serialized    raw bytes (serialized)
##

######################################################################

## 'meta' is a list of functions for updating the file size of the
## database and the file map.

setClass("filehashDB1",
         representation(datafile = "character",
                        meta = "list"),
         contains = "filehash"
         )

setValidity("filehashDB1",
            function(object) {
                    if(!file.exists(object@datafile))
                            return(gettextf("datafile '%s' does not exist",
                                            datafile))
                    TRUE
            })

createDB1 <- function(dbName) {
        if(!hasWorkingFtell())
                stop("need working 'ftell()' to use 'DB1' format")
        if(file.exists(dbName)) {
                message(gettextf("database '%s' already exists", dbName))
                return(TRUE)
        }
        status <- file.create(dbName)

        if(!status)
                stop(gettextf("unable to create database file '%s'", dbName))
        TRUE
}

makeMetaEnv <- function(filename) {
        dbmap <- NULL  ## 'NULL' indicates the map needs to be read
        dbfilesize <- file.info(filename)$size

        updatesize <- function(size) {
                dbfilesize <<- size
        }
        updatemap <- function(map) {
                dbmap <<- map
        }
        getsize <- function() {
                dbfilesize
        }
        getmap <- function() {
                dbmap
        }
        list(updatesize = updatesize,
             updatemap = updatemap,
             getmap = getmap,
             getsize = getsize)
}

initializeDB1 <- function(dbName) {
        if(!hasWorkingFtell())
                stop("need working 'ftell()' to use DB1 format")
        dbName <- normalizePath(dbName)

        new("filehashDB1",
            datafile = dbName,
            meta = makeMetaEnv(dbName),
            name = basename(dbName)
            )
}


readKeyMap <- function(con, map = NULL, pos = 0) {
        if(is.null(map)) {
                ## using 'hash = TRUE' is critical because it can have a major
                ## impact on performance for large databases
                map <- new.env(hash = TRUE, parent = emptyenv())
                pos <- 0
        }
        if(pos < 0)
                stop("'pos' cannot be negative")
        filename <- path.expand(summary(con)$description)
        filesize <- file.info(filename)$size

        if(pos > filesize)
                stop("'pos' cannot be greater than file size")
        .Call("read_key_map", filename, map, filesize, pos)
}

readSingleKey <- function(con, map, key) {
        start <- map[[key]]

        if(is.null(start))
                stop(gettextf("unable to obtain value for key '%s'", key))

        seek(con, start, rw = "read")
        unserialize(con)
}

readKeys <- function(con, map, keys) {
        r <- lapply(keys, function(key) readSingleKey(con, map, key))
        names(r) <- keys
        r
}

gotoEndPos <- function(con) {
        ## Move connection to the end
        seek(con, 0, "end")
        seek(con)
}

writeNullKeyValue <- function(con, key) {
        writestart <- gotoEndPos(con)

        handler <- function(cond) {
                ## Rewind the file back to where writing began and truncate at
                ## that position
                seek(con, writestart, "start", "write")
                truncate(con)
                cond
        }
        tryCatch({
                serialize(key, con)

                len <- as.integer(-1)
                serialize(len, con)
        }, interrupt = handler, error = handler, finally = {
                flush(con)
        })
}

writeKeyValue <- function(con, key, value) {
        writestart <- gotoEndPos(con)

        handler <- function(cond) {
                ## Rewind the file back to where writing began and
                ## truncate at that position; this is probably a bad
                ## idea for files > 2GB
                seek(con, writestart, "start", "write")
                truncate(con)
                cond
        }
        tryCatch({
                serialize(key, con)

                byteData <- serialize(value, NULL)
                len <- length(byteData)
                serialize(len, con)

                writeBin(byteData, con)
        }, interrupt = handler, error = handler, finally = {
                flush(con)
        })
}

setMethod("lockFile", "file", function(db, ...) {
        ## Use 3 underscores for lock file
        sprintf("%s___LOCK", summary(db)$description)
})

createLockFile <- function(name) {
        if(.Platform$OS.type != "windows") 
                status <- .Call("lock_file", name)
        else {
                ## TODO: are these optimal values for max.attempts
                ## and sleep.duration?
                max.attempts <- 4
                sleep.duration <- 0.5
                attempts <- 0
                status <- -1
                while ((attempts <= max.attempts) && ! isTRUE(status >= 0)) {
                        attempts <- attempts + 1
                        status <- .Call("lock_file", name)

                        if(!isTRUE(status >= 0))
                                Sys.sleep(sleep.duration)
                }
        }
        if(!isTRUE(status >= 0))
                stop("cannot create lock file ", sQuote(name))
        TRUE
}

deleteLockFile <- function(name) {
        if(!file.remove(name))
                stop(paste('cannot remove lock file "', name, '"', sep=''))
        TRUE
}

################################################################################
## Internal utilities

filesize <- gotoEndPos

setGeneric("checkMap", function(db, ...) standardGeneric("checkMap"))

setMethod("checkMap", "filehashDB1",
          function(db, filecon, ...) {
                  old.size <- db@meta$getsize()
                  cur.size <- tryCatch({
                          filesize(filecon)
                  }, error = function(err) {
                          old.size
                  })
                  size.change <- old.size != cur.size
                  map <- getMap(db)
                  map0 <- map

                  if(is.null(map))
                          map <- readKeyMap(filecon)
                  else if(size.change) {
                          ## Modify 'map.old' directly
                          map <- tryCatch({
                                  readKeyMap(filecon, map, old.size)
                          }, error = function(err) {
                                  message(conditionMessage(err))
                                  map0
                          })
                  }
                  else
                          map <- map0
                  if(!identical(map, map0)) {
                          db@meta$updatemap(map)
                          db@meta$updatesize(cur.size)
                  }
                  invisible(db)
          })


setGeneric("getMap", function(db) standardGeneric("getMap"))

setMethod("getMap", "filehashDB1",
          function(db) {
                  db@meta$getmap()
          })

################################################################################
## Interface functions

openDBConn <- function(filename, mode) {
        con <- try({
                file(filename, mode)
        }, silent = TRUE)

        if(inherits(con, "try-error"))
                stop("unable to open connection to database")
        con
}

setMethod("dbInsert",
          signature(db = "filehashDB1", key = "character", value = "ANY"),
          function(db, key, value, ...) {
                  con <- openDBConn(db@datafile, "ab")
                  on.exit(close(con))

                  lockname <- lockFile(con)
                  createLockFile(lockname)
                  on.exit(deleteLockFile(lockname), add = TRUE)

                  invisible(writeKeyValue(con, key, value))
          })

setMethod("dbFetch",
          signature(db = "filehashDB1", key = "character"),
          function(db, key, ...) {
                  con <- openDBConn(db@datafile, "rb")
                  on.exit(close(con))

                  lockname <- lockFile(con)
                  createLockFile(lockname)
                  on.exit(deleteLockFile(lockname), add = TRUE)

                  checkMap(db, con)
                  map <- getMap(db)

                  val <- readSingleKey(con, map, key)
                  val
          })

setMethod("dbMultiFetch",
          signature(db = "filehashDB1", key = "character"),
          function(db, key, ...) {
                  con <- openDBConn(db@datafile, "rb")
                  on.exit(close(con))

                  lockname <- lockFile(con)
                  createLockFile(lockname)
                  on.exit(deleteLockFile(lockname), add = TRUE)

                  checkMap(db, con)
                  map <- getMap(db)

                  readKeys(con, map, key)
          })

setMethod("dbExists", signature(db = "filehashDB1", key = "character"),
          function(db, key, ...) {
                  dbkeys <- dbList(db)
                  key %in% dbkeys
          })

setMethod("dbList", "filehashDB1",
          function(db, ...) {
                  con <- openDBConn(db@datafile, "rb")
                  on.exit(close(con))

                  lockname <- lockFile(con)
                  createLockFile(lockname)
                  on.exit(deleteLockFile(lockname), add = TRUE)

                  checkMap(db, con)
                  map <- getMap(db)

                  if(length(map) == 0)
                          character(0)
                  else {
                          keys <- as.list(map, all.names = TRUE)
                          use <- !sapply(keys, is.null)
                          names(keys[use])
                  }
          })

setMethod("dbDelete", signature(db = "filehashDB1", key = "character"),
          function(db, key, ...) {
                  con <- openDBConn(db@datafile, "ab")
                  on.exit(close(con))

                  lockname <- lockFile(con)
                  createLockFile(lockname)
                  on.exit(deleteLockFile(lockname), add = TRUE)

                  invisible(writeNullKeyValue(con, key))
          })

setMethod("dbUnlink", "filehashDB1",
          function(db, ...) {
                  file.remove(db@datafile)
          })

reorganizeDB <- function(db, ...) {
        datafile <- db@datafile

        ## Find a temporary file name
        tempdata <- paste(datafile, "Tmp", sep = "")
        i <- 0
        while(file.exists(tempdata)) {
                i <- i + 1
                tempdata <- paste(datafile, "Tmp", i, sep = "")
        }
        if(!dbCreate(tempdata, type = "DB1")) {
                warning("could not create temporary database")
                return(FALSE)
        }
        on.exit(file.remove(tempdata))

        tempdb <- dbInit(tempdata, type = "DB1")
        keys <- dbList(db)

        ## Copy all keys to temporary database
        nkeys <- length(keys)
        cat("Reorganizing database: ")

        for(i in seq_along(keys)) {
                key <- keys[i]
                msg <- sprintf("%d%% (%d/%d)", round (100 * i / nkeys),
                               i, nkeys)
                cat(msg)

                dbInsert(tempdb, key, dbFetch(db, key))

                back <- paste(rep("\b", nchar(msg)), collapse = "")
                cat(back)
        }
        cat("\n")
        status <- file.rename(tempdata, datafile)

        if(!isTRUE(status)) {
                on.exit()
                warning("temporary database could not be renamed and is left in ",
                        tempdata)
                return(FALSE)
        }
        on.exit()
        cat("Finished; reload database with 'dbInit'\n")
        TRUE
}

setMethod("dbReorganize", "filehashDB1", reorganizeDB)


################################################################################
## Test system's ftell()

hasWorkingFtell <- function() {
        tfile <- tempfile()
        con <- file(tfile, "wb")

        tryCatch({
                bytes <- raw(10)
                begin <- seek(con)

                if(begin != 0)
                        return(FALSE)
                writeBin(bytes, con)
                end <- seek(con)
                offset <- end - begin
                isTRUE(offset == 10)
        }, error = function(e) {
                FALSE
        }, finally = {
                close(con)
                unlink(tfile)
        })
}

######################################################################


