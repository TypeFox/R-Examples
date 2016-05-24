## under Unix-like OS, we don't have choose.files() and Filters
getFileNames <-
function(fPath=".", fNames, fPat) {
    ## do we have file names or a name pattern?
    files <- if(!missing(fNames)) {          # we have file names
        if(is.null(fPath) || is.na(fPath)) { # no path is given
            fNames                           # try file names without path
        } else {                             # we have a path
            paste(path.expand(fPath), fNames, sep="/")
        }
    } else if(!missing(fPat)) {              # we have a name pattern
        list.files(path=fPath, pattern=fPat, full.names=TRUE)
    } else {
        character(0)
    }

    return(files)
}
