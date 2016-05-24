.onAttach <- function(lib, pkg) {
    where <- match(paste("package:", pkg, sep = ""), search())
    ver <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
    ver <- as.character(ver)
    title <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Title")
    title <- as.character(title)
    text <- paste(title, " (version ", ver, ")\n", sep = "")
    packageStartupMessage(text)
}
