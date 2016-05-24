## no S4 methodology here; speedup :
.noGenerics <- TRUE

Sleuth2Manual <- function(){
  viewer <- getOption("pdfviewer")
  file <- system.file("doc/Sleuth2-manual.pdf", package="Sleuth2")
  system(paste(viewer, file, sep=" "))
  invisible()
}

