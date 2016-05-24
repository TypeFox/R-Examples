package.version <- function(all.available=FALSE, lib.loc=NULL)
{
    ##if(is.null(pkg))
    pkg <- .packages(all.available, lib.loc)
    if(is.null(lib.loc)) libpath <- .Library
    pkgs <- NULL
    vers <- NULL
    for(i in pkg) {
        verline <- NULL
        dnm <- paste(libpath, "/", i, sep="")
        flnm <- paste(dnm, "/DESCRIPTION", sep="")
        if(!file.exists(dnm)) {
            cat("\n", i, "not exists.")
            next
        } else if(!file.exists(flnm)) {
            cat("\nDESCRIPTION for", i, "not found.")
            next
        }
        chlines <- readLines(con=flnm, n=-1, ok=TRUE)
        for(j in chlines) 
            verline <- paste(verline, grep("version:", j, ignore.case=T, value=T))
        version <- gsub("Version: ", "", verline, ignore.case=TRUE)
        version <- gsub(" ", "", version)
        cat("\nVersion of", i, ":", version)
        pkgs <- c(pkgs, i)
        vers <- c(vers, version)
    }
    cat("\n")
    invisible(list(package=pkgs, version=vers))
}
