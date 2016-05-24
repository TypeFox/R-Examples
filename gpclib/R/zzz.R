.onAttach <- function(lib, pkg) {
    ver <- read.dcf(file.path(lib,pkg,"DESCRIPTION"), "Version")
    msg <- sprintf("General Polygon Clipper Library for R (version %s)\n\tType 'class ? gpc.poly' for help", as.character(ver))
    packageStartupMessage(msg)
}
