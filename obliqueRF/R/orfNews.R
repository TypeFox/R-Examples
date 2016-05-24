orfNews <- function() {
    newsfile <- file.path(system.file(package="obliqueRF"), "NEWS")
    file.show(newsfile)
}
