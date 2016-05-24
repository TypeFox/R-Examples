##' Convert PDF file to print quality png (default 300 dpi)
##'
##' Access to ghostscript program 'gs' is needed
##' @title Convert pdf to raster format
##' @param files Vector of (pdf-)filenames to process
##' @param dpi DPI
##' @param resolution Resolution of raster image file
##' @param gs Optional ghostscript command
##' @param gsopt Optional ghostscript arguments
##' @param resize Optional resize arguments (mogrify)
##' @param format Raster format (e.g. png, jpg, tif, ...)
##' @param \dots Additional arguments
##' @seealso \code{dev.copy2pdf}, \code{printdev}
##' @export
##' @author Klaus K. Holst
##' @keywords iplot
pdfconvert <- function(files, dpi=300, resolution=1024, gs, gsopt, resize, format="png", ...) {
    if (missing(gsopt))
        gsopt <- "-dSAFTER -dBATCH -dNOPAUSE -sDEVICE=png16m -dGraphicsAlphaBits=4 -dTextAlphaBits=4"
    if (missing(gs)) {
        gs <- names(which(Sys.which(c("gs", "gswin32c", "gswin64c")) != ""))
    }
    cmd1 <- paste0(gs," -r",dpi," -dBackgroundColor='16#ffffff'")
    if (missing(resize)) {
        resize <- paste0("mogrify -resize ", resolution)
    }
    for (f in files) {
        f0 <- strsplit(f,".pdf")[1]
        f.out <- paste(f0,format,sep=".")
        f.pdf <- paste(f0,"pdf",sep=".")
        mycmd1 <- paste0(cmd1, " ", gsopt, " -sOutputFile=", f.out, " > /dev/null ", f.pdf)
        mycmd2 <- paste0(resize, " ", f.out)
        cat(f.pdf)
        system(mycmd1)
        cat(" -> ")
        system(mycmd2)
        cat(f.out, "\n")
  }
}
