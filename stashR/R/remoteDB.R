######################################################################
## Class definitions

## Use 'http://' type URLs
setClass("remoteDB",
         representation(url = "character",
                        dir = "character",
			reposVersion = "numeric"),
	 prototype = list(reposVersion = -1), 
         contains = "filehash"
         )




################################################################################
## Method definitions for 'remoteDB'
################################################################################

setMethod("dbInsert",
          signature(db = "remoteDB", key = "character", value = "ANY"),
          function(db, key, value, ...) {
                  stop("cannot insert into a 'remoteDB' database")
          })

setMethod("dbFetch", signature(db = "remoteDB", key = "character"),
          function(db, key, ...) {
                  if(!dataFileExists(db, key))	
                          getdata(db,key)
                  con <- gzfile(local.file.path(db, key), "rb")
                  on.exit(close(con))
                  unserialize(con) 
          })

setMethod("dbDelete", signature(db = "remoteDB", key = "character"),
          function(db, key, ...) {
                  stop("cannot delete from a 'remoteDB' database")
          })

setMethod("dbList", "remoteDB",
          function(db, ...){
                  info <- reposVersionInfo(db)

                  if(length(info) != 0){
                          keyFiles <- strsplit(info, ":")[[1]][2]
                          keyFilesSep <- strsplit(keyFiles," ", fixed = TRUE)[[1]]

                          ## Strip the object version numbers off
                          gsub("\\.[0-9]+$", "", keyFilesSep)
                  }
                  else
                          character(0)
          })

setMethod("dbExists", signature(db = "remoteDB", key = "character"),
          function(db, key, ...){
                  key %in% dbList(db)  ## returns a vector of TRUE/FALSE
          })


removeCruft <- function(db, keys) {
        keyfiles <- sapply(keys, function(key) local.file.path(db, key))
        sigfiles <- sapply(keys, function(key) local.file.path.SIG(db, key))
        keepFiles <- basename(c(keyfiles, sigfiles))

        fileList <- dir(file.path(db@dir, "data"), full.names = TRUE,
                        all.files = TRUE)
        ## exclude directories
        use <- !file.info(fileList)$isdir
        fileList <- basename(fileList[use])
        removeList <- setdiff(fileList, keepFiles)

        for(filename in removeList) 
                file.remove(file.path(db@dir, "data", filename))
        invisible()
}

cacheVersionFile <- function(db) {  ## 'db' is a 'remoteDB' here
        rpath <- versionFile(db, where = "remote")
        download.file(rpath, versionFile(db, where = "local"), mode = "w",
                      cacheOK = FALSE, quiet = TRUE)
}

setGeneric("dbSync", function(db, ...) standardGeneric("dbSync"))

setMethod("dbSync", signature(db = "remoteDB"),
          function(db, key = NULL, ...) {
                  if(db@reposVersion != -1)
                          stop("no synchronization for a 'remoteDB' object ",
                               "with a fixed version")
                  cacheVersionFile(db)
                  
                  if(!is.null(key)) {
                          isLocal <- dataFileExists(db, key)

                          if(!all(isLocal))
                                  stop("not all files referenced in the 'key' ",
                                       "vector were previously downloaded, ",
                                       "no files updated")
                  }
                  else {
                          key <- dbList(db)
                          isLocal <- dataFileExists(db, key)
                  }
                  use <- which(isLocal)
                  key <- names(isLocal[use])

                  ## Get the latest versions of all local keys
                  for(i in seq(along = key)) 
                          dbFetch(db, key[i])

                  ## Remove cruft from previous versions
                  removeCruft(db, key)
          })




########################################################################
## Utility Functions  ##################################################
########################################################################

setGeneric("versionFile", function(db, ...) standardGeneric("versionFile"))

## Return the URL for the 'version' file
setMethod("versionFile", "remoteDB",
          function(db, where = c("remote", "local"), ...) {
                  where <- match.arg(where)

                  switch(where,
                         remote = paste(db@url, "version", sep = "/"),
                         local = file.path(db@dir, "version")
                         )
          })


## use readLines to find last line or the line corresponding to
## 'db@reposVersion', returns as character string

readVersionFileLine <- function(db) {
        ## Always read the local copy of the version file.  The number of
        ## lines read is equal to the version of the repository.  If the
        ## repository version is -1, then we read the whole file
        vfile <- versionFile(db, where = "local")

        tryCatch({
                VerList <- readLines(vfile, n = db@reposVersion)
                VerList[length(VerList)]
        }, error = function(cond) {
                character(0)
        })
}

setGeneric("reposVersionInfo",
           function(db, ...) standardGeneric("reposVersionInfo"))

setMethod("reposVersionInfo", "remoteDB",
          function(db, ...) {
                  readVersionFileLine(db)
          })

## returns "object version" associated with a given key
## 
## Get the version number for an object corresponding to
## 'db@reposVersion'

getSpecificObjectVersion <- function(db, keyList) {
        info <- reposVersionInfo(db)
        currNum <- numeric(length(keyList))
        
        if(length(info) == 0) 
                return(currNum)
        keyFiles <- strsplit(info, ":")[[1]][2]
        keyFilesSep <- strsplit(keyFiles, " ", fixed = TRUE)[[1]]

        sapply(keyList, function(key) {
                v <- grep(paste("^", key, "\\.[0-9]+$", sep = ""),
                          keyFilesSep, value = TRUE)
                if(length(v) > 0)
                        as.numeric(substring(v, nchar(key) + 2))
                else
                        0
        }, USE.NAMES = FALSE)
}

setGeneric("objectVersion", function(db, ...) standardGeneric("objectVersion"))

setMethod("objectVersion", "remoteDB",
          function(db, key, ...) {
                  getSpecificObjectVersion(db, key)
          })




## function returning the integer corresponding to the lastest repos
## version by reading first integer of the last line of the version
## file

latestReposVersionNum <- function(info){ 
        if(length(info) != 0) 
                as.numeric(strsplit(info, ":")[[1]][1])
        else
                0
}


## reads last line of the version file creates new info on repository
## version to be inserted in last line of the version file updates
## keepKey = FALSE deletes the key from the version information


updatedReposVersionInfo <- function(db, key, keepKey = TRUE){
        ## find and remove key (for updating) from latest repository version info ##
        info <- reposVersionInfo(db)
        reposV <- latestReposVersionNum(info) + 1

        if(length(info) != 0){
                keyFiles <- strsplit(info, ":")[[1]][2]
                keyFilesSep <- strsplit(keyFiles," ")[[1]]
                v <- grep(paste("^",key,"\\.[0-9]+$",sep=""),keyFilesSep)

                if(length(v)==0)
                        others <- keyFilesSep
                else
                        others <- keyFilesSep[-v]
                if(length(others)==0)
                        othersCollapsed <- paste(others, collapse=" ")
                else{	
                        if(length(others)==1 & is.na(others[1]))
                                othersCollapsed <- character(0)
                        else 
                                othersCollapsed <- paste(others, collapse=" ")
                }
                newKeyInfo <- paste(key,objectVersion(db,key)+1, sep=".")
                if(keepKey)
                        updatedKeyFiles <- paste(othersCollapsed, newKeyInfo,
                                                 sep = " ")
                else
                        updatedKeyFiles <- othersCollapsed
        }
        else
                updatedKeyFiles <- paste(key,1,sep=".")
        ## remove leading space resulting from intially inserting same key twice#
        updatedKeyFiles <- gsub("^ ","",updatedKeyFiles)
        paste(reposV, updatedKeyFiles, sep = ":") 
}



######### updating version file ########## 
updateVersion <- function(db, key, keepKey = TRUE){
        cat(updatedReposVersionInfo(db,key,keepKey),
            file = versionFile(db), sep = "\n", append = TRUE)
}


###############################
## local.file.path ############  Creates a file path in the local data  
###############################  directory (to be used internally).	

key2filename <- function(key) {
        if(!is.character(key))
                stop("'key' should be character")
        sapply(key, digest, serialize = FALSE, USE.NAMES = FALSE)
}

local.file.path <- function(db, key, objVerNum = objectVersion(db, key)) {
        key.filename <- paste(key2filename(key), objVerNum, sep = ".")
        file.path(db@dir, "data", key.filename)
}

###############################
## local.file.path.SIG ######## Creates a file path in the local data  
############################### directory (to be used internally) for the SIG files.

local.file.path.SIG <- function(db, key, objVerNum = objectVersion(db, key)) {
        key.SIGfilename <- paste(key2filename(key), objVerNum, "SIG", sep = ".")
        file.path(db@dir, "data", key.SIGfilename)
}



#################### 
## dataFileExists ##
#################### 

## Returns TRUE if data file for 'key' is in local dir, otherwise
## returns FALSE. We have 'key' allowed to be a character vector with
## more than one key.

dataFileExists <- function(db, key) {
        datadir <- file.path(db@dir, "data")

        if(!file.exists(datadir))
                stop("local data directory does not exist")

        ver <- objectVersion(db, key)
        val <- file.exists(local.file.path(db, key, ver))
        names(val) <- key
        val
}

####################
## getdata ######### downloads the key & the SIG file
####################

remote.file.path <- function(db, key) {
        key.filename <- paste(key2filename(key), objectVersion(db, key),
                              sep = ".")
        file.path(db@url, "data", key.filename)
}

remote.file.path.SIG <- function(db, key) {
        key.SIGfilename <- paste(key2filename(key), objectVersion(db, key), "SIG",
                                 sep = ".")
        file.path(db@url, "data", key.SIGfilename)
}

getdata <- function(db, key) {
        localFiles <- c(data = local.file.path(db, key),
                        sig = local.file.path.SIG(db, key))

        handler <- function(cond) {
                ## If a condition is thrown (e.g. error or interrupt), delete
                ## whatever was downloaded
                file.remove(localFiles)
                cond
        }
        status <- tryCatch({
                remotePath <- remote.file.path(db, key)
                download.file(remotePath, localFiles["data"], mode = "wb",
                              cacheOK = FALSE, quiet = stashROption("quietDownload"))

                remoteSIGPath <- remote.file.path.SIG(db, key)
                download.file(remoteSIGPath, localFiles["sig"], mode = "wb",
                              cacheOK = FALSE, quiet = stashROption("quietDownload"))
        }, error = handler, interrupt = handler)

        if(inherits(status, "condition")) 
                stop(gettextf("problem downloading data for key '%s'", key))
        invisible(status)
}




