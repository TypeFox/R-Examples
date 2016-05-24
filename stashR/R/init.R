## 'dbName' here should be a character vector with names 'url' and
## optionally 'dir', e.g. c(url = "http://....", dir = "/home/rpeng").

setMethod("initialize", "remoteDB",
          function(.Object, ...) {
                  .Object <- callNextMethod()
                  .Object <- dbCreate(.Object)
                  cacheVersionFile(.Object)
                  .Object
          })

setMethod("initialize", "localDB",
          function(.Object, ...) {
                  .Object <- callNextMethod()
                  .Object <- dbCreate(.Object)
                  .Object
          })

createLocalDir <- function(db) {
        ## create the local main directory and data sub-directory to
        ## store the data files ##
        datadir <- file.path(db@dir,"data")
        status <- dir.create(db@dir, showWarnings = FALSE,
                             recursive = TRUE)
        
        if(!status && !file.exists(db@dir))
                stop(gettextf("problem creating directory '%s'", db@dir))
        status <- dir.create(datadir, showWarnings = FALSE,
                             recursive = TRUE)
        
        if(!status && !file.exists(datadir))
                stop(gettextf("problem creating directory '%s'", datadir))
}

setMethod("dbCreate",
          signature(db = "remoteDB"),
          function(db, ...) {
                  ## remove trailing "/" on dir and url ##
                  db@dir <- sub("\\/$","", db@dir, perl = TRUE)
                  db@url <- sub("\\/$","", db@url, perl = TRUE)
                  
                  createLocalDir(db)
                  
                  ## save url in the R workspace format in the main directory ##
                  myurl <- db@url 
                  save(myurl, file = file.path(db@dir, "url"))
                  invisible(db)
          })

setMethod("dbCreate",
          signature(db = "localDB"),
          function(db, ...) {
                  ## remove trailing "/" on dir ##
                  db@dir <- sub("\\/$","", db@dir, perl = TRUE)

                  createLocalDir(db)

                  if(!file.exists(versionFile(db)))
                          file.create(versionFile(db))
                  invisible(db)
          })




