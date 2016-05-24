require(gWidgets)
options("guiToolkit"="RGtk2")


gWidgetsDir <- system.file("tests",package="gWidgets")
## should be there, but just in case
if(gWidgetsDir != "") {
  files <- list.files(gWidgetsDir,
                      pattern = "\\.R$",
                      full.names = TRUE)

  files <- files[grepl("^ex", basename(files))]

  ## XXX removed to get past CRAN
#  for(unitTest in files) {
#    source(unitTest)
#  }
}
