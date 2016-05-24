.onAttach <- function(lib, pkg) {
    dcf <- read.dcf(file.path(lib, pkg, "DESCRIPTION"))
    msg <- gettextf("%s (%s %s)", dcf[, "Title"],
                    as.character(dcf[, "Version"]), dcf[, "Date"])
    packageStartupMessage(paste(strwrap(msg), collapse = "\n"))
}

.onLoad <- function(lib, pkg) {
    if(!require(boot))
        stop("'boot' package required")
}

