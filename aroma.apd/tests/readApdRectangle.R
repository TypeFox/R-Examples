library("aroma.apd")

# Local functions
rotate270 <- function(x, ...) {
  x <- t(x)
  nc <- ncol(x)
  if (nc < 2) return(x)
  x[,nc:1,drop=FALSE]
}

# Scan current directory for APD files
files <- list.files(pattern="[.](apd|APD)$")
files <- files[!file.info(files)$isdir]
if (length(files) > 0) {
  apdFile <- files[1]

  # Read APD intensities in the upper left corner
  apd <- readApdRectangle(apdFile, xrange=c(0,250), yrange=c(0,250))
  z <- rotate270(apd$intensities)
  sub <- paste("Chip type:", apd$header$chipType)
  image(z, col=gray.colors(256), axes=FALSE, main=apdFile, sub=sub)
  text(x=0, y=1, labels="(0,0)", adj=c(0,-0.7), cex=0.8, xpd=TRUE)
  text(x=1, y=0, labels="(250,250)", adj=c(1,1.2), cex=0.8, xpd=TRUE)
}
