"lmomcoBook" <-
function() {
    file <- file.path(system.file(package="lmomco"), "ERRATA_FOR_ISBN9781463508418.txt")
    file.show(file)
}

"lmomcoNews" <-
function() {
    file <- file.path(system.file(package="lmomco"), "NEWS")
    file.show(file)
}
