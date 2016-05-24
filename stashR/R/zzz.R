.onAttach <- function(lib, pkg) {
        dcf <- read.dcf(file.path(lib, pkg, "DESCRIPTION"))
        msg <- gettextf("%s (%s %s)", dcf[, "Title"],
                        as.character(dcf[, "Version"]), dcf[, "Date"])
        packageStartupMessage(paste(strwrap(msg), collapse = "\n"))
}

.onLoad <- function(lib, pkg) {
        if(!capabilities("http/ftp"))
                warning("'http/ftp' capabilities not available")
        stashROption("quietDownload", FALSE)
}

.stashROptions <- new.env()

## Valid options:
##
## quietDownload:  Should download progress be shown?

stashROption <- function(name, value) {
        if(missing(name))
                as.list(.stashROptions)
        else if(missing(value))
                get(name, .stashROptions, inherits = FALSE)
        else
                assign(name, value, .stashROptions, inherits = FALSE)
}
