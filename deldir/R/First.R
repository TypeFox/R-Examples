.onLoad <- function(lib, pkg) {
	library.dynam("deldir", pkg, lib)
}

.onAttach <- function(lib, pkg) {
	ver <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
        packageStartupMessage(paste(pkg, ver))
}
