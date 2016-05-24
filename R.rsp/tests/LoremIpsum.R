library("R.rsp")

path <- system.file(package="R.rsp")
path <- file.path(path, "rsp_LoremIpsum")
pathnames <- list.files(path=path, pattern="[.]rsp$", full.names=TRUE)
for (pathname in pathnames) {
  outPath <- gsub("LoremIpsum.", "", basename(pathname))
  outPath <- file.path("LoremIpsum", outPath)
  tryCatch({
    pathnameR <- rfile(pathname, workdir=outPath, verbose=-10)
    print(pathnameR)
  }, error = function(ex) {
    print(ex);
  })
}
