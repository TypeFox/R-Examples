.dbEnv <- new.env()

.onLoad <- function(lib, pkg) {
    stashROption("quietDownload", TRUE)
}

.onAttach <- function(lib, pkg) {
    dcf <- read.dcf(file.path(lib, pkg, "DESCRIPTION"))
    msg <- gettextf("%s (%s %s)", dcf[, "Title"],
                    as.character(dcf[, "Version"]), dcf[, "Date"])
    packageStartupMessage(paste(strwrap(msg), collapse = "\n"))
    packageStartupMessage("Initialize database using 'initMCAPS'")
}

