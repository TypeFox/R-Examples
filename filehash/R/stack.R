setClass("stack",
         representation(stack = "filehashDB1",
                        name = "character"))

setMethod("show", "stack",
          function(object) {
                  cat(gettextf("<stack: %s>\n", object@name))
                  invisible(object)
          })

createS <- function(filename) {
        dbCreate(filename, "DB1")
        stack <- dbInit(filename, "DB1")
        dbInsert(stack, "top", NULL)

        new("stack", stack = stack, name = filename)
}

initS <- function(filename) {
        new("stack",
            stack = dbInit(filename, "DB1"),
             name = filename)
}

setMethod("lockFile", "stack",
          function(db, ...) {
                  paste(db@name, "slock", sep = ".")
          })

setMethod("push", c("stack", "ANY"), function(db, val, ...) {
        node <- list(value = val,
                    nextkey = dbFetch(db@stack, "top"))
        topkey <- sha1(node)

        createLockFile(lockFile(db))
        on.exit(deleteLockFile(lockFile(db)))

        dbInsert(db@stack, topkey, node)
        dbInsert(db@stack, "top", topkey)
})

setGeneric("mpush", function(db, vals, ...) standardGeneric("mpush"))

setMethod("mpush", c("stack", "ANY"), function(db, vals, ...) {
        if(!is.list(vals))
                vals <- as.list(vals)
        createLockFile(lockFile(db))
        on.exit(deleteLockFile(lockFile(db)))

        topkey <- dbFetch(db@stack, "top")

        for(i in seq_along(vals)) {
                node <- list(value = vals[[i]],
                             nextkey = topkey)
                topkey <- sha1(node)

                dbInsert(db@stack, topkey, node)
                dbInsert(db@stack, "top", topkey)
        }
})

setMethod("isEmpty", "stack", function(db, ...) {
        h <- dbFetch(db@stack, "top")
        is.null(h)
})


setMethod("top", "stack", function(db, ...) {
        createLockFile(lockFile(db))
        on.exit(deleteLockFile(lockFile(db)))

        if(isEmpty(db))
                stop("stack is empty")
        h <- dbFetch(db@stack, "top")
        node <- dbFetch(db@stack, h)
        node$value
})

setMethod("pop", "stack", function(db, ...) {
        createLockFile(lockFile(db))
        on.exit(deleteLockFile(lockFile(db)))

        if(isEmpty(db))
                stop("stack is empty")
        h <- dbFetch(db@stack, "top")
        node <- dbFetch(db@stack, h)

        dbInsert(db@stack, "top", node$nextkey)
        dbDelete(db@stack, h)
        node$value
})
