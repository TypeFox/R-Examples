aCRMNews <-
function() {
    newsfile <- file.path(system.file(package="aCRM"), "NEWS")
    file.show(newsfile)
}
