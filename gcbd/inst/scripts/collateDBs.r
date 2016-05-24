#!/usr/bin/r -t

suppressMessages(library(RSQLite))

readFirst <- function(dbfile) {
    if (file.exists(dbfile)) {
        dbcon <- dbConnect(dbDriver("SQLite"), dbname=dbfile)
        d <- dbReadTable(dbcon, "benchmark")
        dbDisconnect(dbcon)
    }
    invisible(d)
}

appendToSecond <- function(newdf, dbfile) {
    if (file.exists(dbfile)) {
        dbcon <- dbConnect(dbDriver("SQLite"), dbname=dbfile)
        res <- dbWriteTable(dbcon, "benchmark", newdf, row.names=FALSE, overwrite=FALSE, append=TRUE)
        dbDisconnect(dbcon)
    }
    invisible(NULL)
}


if (length(argv) != 2) stop("Need two arguments!\n")

d <- readFirst(argv[1])
appendToSecond(d, argv[2])
