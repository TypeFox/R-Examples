.packageName <- "retimes"

.onAttach <- function(lib, pkg) {
    where <- match(paste("package:", pkg, sep = ""), search())
    ver <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
    ver <- as.character(ver)
    title <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Title")
    title <- as.character(title)
    packageStartupMessage(paste(title, " (version ", ver, ")", sep = ""))
    #library.dynam("retimes", pkg, lib)
}

.onLoad <- function(lib, pkg)
    library.dynam("retimes", pkg, lib)
