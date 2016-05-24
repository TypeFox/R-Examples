#!/usr/bin/r -t

suppressMessages(library(RSQLite))

runCrosstabs <- function(dbfile) {
    if (file.exists(dbfile)) {
        dbcon <- dbConnect(dbDriver("SQLite"), dbname=dbfile)
        d <- dbReadTable(dbcon, "benchmark")
        dbDisconnect(dbcon)

        cat("\nFor", dbfile, "\n")
        print(xtabs( ~ type + host, data=d))
    }
    invisible(NULL)
}

runCrosstabs("/var/tmp/gcbd.sqlite")
runCrosstabs("../sql/gcbd.sqlite")
