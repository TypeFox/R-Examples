## no S4 methodology here; speedup :
.noGenerics <- TRUE

Sleuth3Manual <- function(){
  viewer <- getOption("pdfviewer")
  file <- system.file("doc/Sleuth3-manual.pdf", package="Sleuth3")
  system(paste(viewer, file, sep=" "))
  invisible()
}

