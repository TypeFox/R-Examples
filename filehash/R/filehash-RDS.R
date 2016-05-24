######################################################################
## Copyright (C) 2006, Roger D. Peng <rpeng@jhsph.edu>
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

################################################################################
## Class 'filehashRDS'

setClass("filehashRDS",
         representation(dir = "character"),
         contains = "filehash"
         )

setValidity("filehashRDS",
            function(object) {
                    if(length(object@dir) != 1)
                            return("only one directory should be set in 'dir'")
                    if(!file.exists(object@dir))
                            return(gettextf("directory '%s' does not exist",
                                            object@dir))
                    TRUE
            })

createRDS <- function(dbName) {
        if(!file.exists(dbName)) {
                status <- dir.create(dbName)

                if(!status)
                        stop(gettextf("unable to create database directory '%s'",
                                      dbName))
        }
        else
                message(gettextf("database '%s' already exists", dbName))
        TRUE
}

initializeRDS <- function(dbName) {
        ## Trailing '/' causes a problem in Windows?
        dbName <- sub("/$", "", dbName, perl = TRUE)
        new("filehashRDS", dir = normalizePath(dbName),
            name = basename(dbName))
}

## For case-insensitive file systems, objects with the same name but
## differ by capitalization might get clobbered.  `mangleName()'
## inserts a "@" before each capital letter and `unMangleName()'
## reverses the operation.

mangleName <- function(oname) {
        if(any(grep("@",oname,fixed=TRUE))) 
                stop("RDS format cannot cope with objects with @ characters",
                        " in their names")
        gsub("([A-Z])", "@\\1", oname, perl = TRUE)
}

unMangleName <- function(mname) {
        gsub("@", "", mname, fixed = TRUE)
}

## Function for mapping a key to a path on the filesystem
setGeneric("objectFile", function(db, key) standardGeneric("objectFile"))
setMethod("objectFile", signature(db = "filehashRDS", key = "character"),
          function(db, key) {
                  file.path(db@dir, mangleName(key))
          })

################################################################################
## Interface functions

setMethod("dbInsert",
          signature(db = "filehashRDS", key = "character", value = "ANY"),
          function(db, key, value, safe = TRUE, ...) {
                  writefile <- if(safe)
                          tempfile()
                  else
                          objectFile(db, key)
                  con <- gzfile(writefile, "wb")

                  writestatus <- tryCatch({
                          serialize(value, con)
                  }, condition = function(cond) {
                          cond
                  }, finally = {
                          close(con)
                  })
                  if(inherits(writestatus, "condition"))
                          stop(gettextf("unable to write object '%s'", key))
                  if(!safe)
                          return(invisible(!inherits(writestatus, "condition")))

                  cpstatus <- file.copy(writefile, objectFile(db, key),
                                        overwrite = TRUE)

                  if(!cpstatus)
                          stop(gettextf("unable to insert object '%s'", key))
                  else {
                          rmstatus <- file.remove(writefile)

                          if(!rmstatus)
                                  warning("unable to remove temporary file")
                  }
                  invisible(cpstatus)
          })

setMethod("dbFetch", signature(db = "filehashRDS", key = "character"),
          function(db, key, ...) {
                  ## Create filename from key
                  ofile <- objectFile(db, key)
                  ## Open connection
                  val <- tryCatch({
                          con<-gzfile(ofile)
                          # note it is necessary to split creating and opening
                          # the connection into two steps so that the connection
                          # can be closed/destroyed successfully if ofile does 
                          # not exist (avoiding connection leaks).
                          open(con,"rb")
                          ## Read data
                          unserialize(con)
                  }, condition = function(cond) {
                          cond
                  }, finally = {
                          close(con)
                  })
                  if(inherits(val, "condition")) 
                          stop(gettextf("unable to obtain value for key '%s'",
                                        key))
                  val
          })

setMethod("dbMultiFetch",
          signature(db = "filehashRDS", key = "character"),
          function(db, key, ...) {
                  r <- lapply(key, function(k) dbFetch(db, k))
                  names(r) <- key
                  r
          })

setMethod("dbExists", signature(db = "filehashRDS", key = "character"),
          function(db, key, ...) {
                  key %in% dbList(db)
          })

setMethod("dbList", "filehashRDS",
          function(db, ...) {
                  ## list all keys/files in the database
                  fileList <- dir(db@dir, all.files = TRUE, full.names = TRUE)
                  use <- !file.info(fileList)$isdir
                  fileList <- basename(fileList[use])

                  unMangleName(fileList)
          })

setMethod("dbDelete", signature(db = "filehashRDS", key = "character"),
          function(db, key, ...) {
                  ofile <- objectFile(db, key)

                  ## remove/delete the file
                  status <- file.remove(ofile)
                  invisible(isTRUE(all(status)))
          })

setMethod("dbUnlink", "filehashRDS",
          function(db, ...) {
                  ## delete the entire database directory
                  d <- db@dir
                  status <- unlink(d, recursive = TRUE)
                  invisible(status)
          })

