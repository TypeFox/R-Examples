WebPlot <- function(width=640, height=480, type='png', ...) {
  file <- paste('tmp-',paste(sprintf('%x',as.integer(runif(4)*65536)),collapse=''),'.tmp',sep='')
  Cairo(width, height, type=type, file=file, ...)
  mime <- switch(type, png="image/png", pdf="application/pdf", jpg="image/jpeg", jpeg="image/jpeg", gif="image/gif", "application/octet")
  structure(list(file=file,type=type,mime=mime,width=width,height=height), class="WebPlot")
}

as.WebResult.WebPlot <- function(x, ...) {
  dev.off()
  structure(c("tmpfile", x$file, x$mime), class="WebResult")
}

as.character.WebPlot <- function(x, ...) {
  dev.off()
  sz <- file.info(x$file)$size
  r <- readBin(x$file, raw(), sz)
  unlink(x$file)
  paste0("<img src='",base64enc::dataURI(r, x$mime),"' width=",x$width," height=",x$height,">")
}
