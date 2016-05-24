rrfNews <- function() {
    newsfile <- file.path(system.file(package="RRF"), "NEWS")
    file.show(newsfile)
}
