snpRFNews <- function() {
    newsfile <- file.path(system.file(package="snpRF"), "NEWS")
    file.show(newsfile)
}
