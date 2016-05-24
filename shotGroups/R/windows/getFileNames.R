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
    } else if(interactive()) {
        ## we are under Windows since file sits in a platform-specific directory
        ## -> we have choose.files() to choose files interactively
        myFilt <- rbind(Filters, txtCsvDat=c("Data files (*.txt, *.csv, *.dat)",
                                             "*.txt;*.csv;*.dat"))
        choose.files(filters=myFilt[c("txtCsvDat", "All"), ], index=1)
    } else {
        character(0)
    }

    return(files)
}

## are we are in interactive mode AND under Windows?
##    } else if(interactive() && (.Platform$OS.type == "windows")) {
