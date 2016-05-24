AUCNews <-
function() {
    newsfile <- file.path(system.file(package="AUC"), "NEWS")
    file.show(newsfile)
}
