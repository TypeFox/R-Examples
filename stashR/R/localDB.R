################################################################################
## 'localDB' classes and methods
################################################################################

## For local directories
setClass("localDB",
         representation(dir = "character",
			reposVersion = "numeric"),
	 prototype = list(reposVersion = -1),
         contains = "filehash"
         )


################################################################################
## Methods for 'localDB'
################################################################################

setMethod("dbUnlink", signature(db = "localDB"),
          function(db, ...) {
                  unlink(db@dir, recursive = TRUE)
          })

## new dbInsert ## (note: I removed overwrite option)

setMethod("dbInsert",
          signature(db = "localDB", key = "character", value = "ANY"),
          function(db, key, value, ...) {
                  if(db@reposVersion != -1) 
                          stop("inserting key into pervious version not allowed")
                  
                  ## update the data files ##
                  vn <- objectVersion(db, key) + 1
                  
                  ## update the 'version' file ##
                  updateVersion(db,key)

                  ## write out 'value' to a file
                  con <- gzfile(local.file.path(db, key, vn), "wb")
                  tryCatch({
                          serialize(value, con)
                  }, finally = {
                          if(isOpen(con))
                                  close(con)
                  })
                  ## 'con' must be closed before taking the digest
                  
                  digest <- md5sum(local.file.path(db, key, vn))
                  digestWithKey <- paste(as.character(digest), key, vn,
                                         sep = "  ")
                  writeLines(digestWithKey,
                             con = local.file.path.SIG(db, key, vn))
          })


setMethod("dbFetch", signature(db = "localDB", key = "character"),
          function(db, key, ...) {
                  if(!dataFileExists(db, key)) 
                          stop(gettextf("key '%s' not in database", key))
                  con <- gzfile(local.file.path(db, key), "rb")
                  on.exit(close(con))
                  unserialize(con)
          })

## doesn't delete files from repository, just deletes key from latest line
## of the version file
setMethod("dbDelete", signature(db = "localDB", key = "character"),
          function(db, key, ...) {
                  if(db@reposVersion != -1) 
                          stop("deleting key from previous version not allowed")
                  if(!(key %in% dbList(db)))
                          stop(gettextf("key '%s' not in current version", key))
                  
                  updateVersion(db, key, keepKey = FALSE)
          })


setMethod("dbList", "localDB",
          function(db, ...){
                  info <- reposVersionInfo(db)

                  if(length(info)!=0){
                          keyFiles <- strsplit(info, ":")[[1]][2]

                          if(!is.na(keyFiles)){
                                  keyFilesSep <- strsplit(keyFiles," ")[[1]]
                                  gsub("\\.[0-9]+$", "", keyFilesSep)
                          }
                          else
                                  character(0)
                  }
                  else
                          character(0)
          })

setMethod("dbExists", signature(db = "localDB", key = "character"),
          function(db, key, ...){
                  key %in% dbList(db)	# returns a vector of T/F
          })

## db <- new("localDB",dir="testlocal",name="test")

## db <- new("remoteDB",url="http://www.biostat.jhsph.edu/~seckel/remoteDBExampleVersioning", dir="remotelocal",name="remote")


################################################################################
## Utilities for 'localDB' objects
################################################################################

## 'getKeyFiles' uses the 'version' file instead of reading the 'data'
## directory directly
## 
## for localDB: 
## determine last version of the object in the repository    ##
## (note this object may have been previously deleted, so    ##
## we need to look in the data directory at the data files)  ##
##
## for remoteDB:
## read pertinent line of the version file from the internet ##

sortByVersionNumber <- function(keyFiles) {
        splitFiles <- strsplit(keyFiles, ".", fixed = TRUE)

        ## The object version number is always at the end
        num <- sapply(splitFiles, function(x) x[length(x)])
        num <- as.numeric(num)
        keyFiles[order(num, decreasing = FALSE)]
}

getKeyFiles <- function(db) {
        version <- readLines(versionFile(db))

        ## Strip repository version numbers
        v <- sub("^[0-9]+:", "", version)
        keyFiles <- unlist(strsplit(v, " ", fixed = TRUE), use.names = FALSE)

        if(is.null(keyFiles))
                keyFiles <- character(0)
        keyFiles
}

## For 'db@reposVersion == -1' in a 'localDB', figure out the latest
## version number for an object

calculateLatestObjectVersion <- function(db, keyList) {
        keyFiles <- getKeyFiles(db)

        sapply(keyList, function(key) {
                use <- grep(paste("^", key, "\\.[0-9]+$", sep = ""), keyFiles)
                oFiles <- sortByVersionNumber(keyFiles[use])
                latestFile <- oFiles[length(oFiles)]
                
                if(length(latestFile) != 0){
                        latestFileSplit <- strsplit(latestFile,".", fixed=TRUE)[[1]]
                        lastObjVer <- latestFileSplit[length(latestFileSplit)]
                        as.numeric(lastObjVer)
                }
                else
                        0
        }, USE.NAMES = FALSE)
}


setMethod("objectVersion", "localDB",
          function(db, key, ...) {
                  if(db@reposVersion == -1)
                          calculateLatestObjectVersion(db, key)
                  else 
                          getSpecificObjectVersion(db, key)
          })


## Return the path for the 'version' file
setMethod("versionFile", "localDB",
          function(db, where = c("local", "remote"), ...) {
                  where <- match.arg(where)
                  switch(where,
                         local = file.path(db@dir, "version"),
                         remote = {
                                 stop("'localDB' databases do not have ",
                                      "remote version files")
                         })
          })


setMethod("reposVersionInfo", "localDB",
          function(db, ...) {
                  readVersionFileLine(db)
          })

