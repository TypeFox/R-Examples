LPSNews <- function() {
    newsfile <- file.path(system.file(package="LPStimeSeries"), "NEWS")
    file.show(newsfile)
}
