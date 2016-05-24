baseurl <- "http://www.biostat.jhsph.edu/MCAPS/v0.3"


initMCAPS <- function(basedir = "MCAPS") {
    if(!file.exists(basedir))
        message(gettextf("creating directory '%s' for local storage", basedir))
    MCAPSinfo <- new("remoteDB",
                     url = paste(baseurl, "MCAPSinfo", sep = "/"),
                     dir = file.path(basedir, "MCAPSinfo"),
                     name = "MCAPSinfo")
    assign("MCAPSinfo", MCAPSinfo, .dbEnv)
}

getData <- function(name = NULL) {
    checkEnv()

    if(is.null(name)) {
        for(dbname in ls(.dbEnv)) {
            objects <- dbList(get(dbname, .dbEnv))
            msg <- paste(dbname, ": ", paste(objects, collapse = ", "),
                         "\n", sep = "")
            writeLines(strwrap(msg))
        }
        invisible(objects)
    }
    else
        dbFetch(.dbEnv$MCAPSinfo, name)
}

checkEnv <- function() {
    if(!exists(".dbEnv") || length(ls(.dbEnv)) == 0)
        stop("need to call 'initMCAPS' first")
    TRUE
}
