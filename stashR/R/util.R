setGeneric("setDir<-", function(db, value) standardGeneric("setDir<-"))

setReplaceMethod("setDir", signature(db = "remoteDB", value = "character"),
                 function(db, value) {
                         db@dir <- value
                         db
                 })

setGeneric("copyDB", function(db, ...) standardGeneric("copyDB"))

setMethod("copyDB", "localDB",
          function(db, dir, ...) {
                  if(identical(dir, db@dir))
                          stop("cannot copy 'localDB' to original directory")
                  newdb <- new("localDB", dir = dir, name = db@name)
                  keys <- dbList(db)

                  for(key in keys) 
                          dbInsert(newdb, key, dbFetch(db, key))
                  newdb
          })

getCurrentReposVersion <- function(db) {
        info <- reposVersionInfo(db)

        if(length(info) > 0) {
                num <- strsplit(info, ":", fixed = TRUE)[[1]][1]
                as.numeric(num)
        }
        else
                0
}

setGeneric("reposVersion",
           function(db, ...) standardGeneric("reposVersion"))

setMethod("reposVersion", "localDB",
          function(db, ...) {
                  getCurrentReposVersion(db)
          })

setMethod("reposVersion", "remoteDB",
          function(db, ...) {
                  getCurrentReposVersion(db)
          })

setGeneric("reposVersion<-",
           function(db, value) standardGeneric("reposVersion<-"))

setReplaceMethod("reposVersion",
                 signature(db = "remoteDB", value = "numeric"),
                 function(db, value) {
                         db@reposVersion <- value
                         db
                 })

setReplaceMethod("reposVersion",
                 signature(db = "localDB", value = "numeric"),
                 function(db, value) {
                         db@reposVersion <- value
                         db
                 })

convertOldStashR <- function() {
        infiles <- dir(".", all.files = TRUE)
        use <- !file.info(infiles)$isdir
        infiles <- infiles[use]
        outfiles <- sub("\\.SIG$", ".1.SIG", infiles)
        sigfiles <- grep("\\.SIG$", outfiles)
        outfiles[-sigfiles] <- paste(outfiles[-sigfiles], "1", sep = ".")

        for(i in seq(along = infiles)) {
                file.rename(infiles[i], outfiles[i])
        }
        keys <- infiles[-sigfiles]
        vline <- paste("1", paste(keys, "1", sep = ".", collapse = " "), sep = ":")
        writeLines(vline, "../version")

}


################################################################################
## Check validity of the data using MD5 digests
################################################################################

readRemoteSIG <- function(db, key) {
        SIGfile <- basename(local.file.path.SIG(db, key))
        con <- url(file.path(db@url, "data", SIGfile))
        open(con, "r")  ## SIG files are text
        on.exit(close(con))

        val <- scan(con, quiet = TRUE, what = "character", sep = " ")[1]
        as.character(val)
}

readLocalSIG <- function(db, key) {
        path <- local.file.path.SIG(db, key)
        val <- scan(path, quiet = TRUE, what = "character", sep = " ")[1]
        as.character(val)
}

validDataRemote <- function(db, key) {
        localSIG <- readLocalSIG(db, key)
        remoteSIG <- readRemoteSIG(db, key)

        isTRUE(localSIG == remoteSIG)
}        

validDataInternal <- function(db, key) {
        localSIG <- readLocalSIG(db, key)
        digest <- md5sum(local.file.path(db, key))
        digest <- as.character(digest)
        isTRUE(digest == localSIG)
}
