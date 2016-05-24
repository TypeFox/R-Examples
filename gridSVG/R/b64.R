base64enc <- function(filename) {
    filesize <- file.info(filename)$size

    if (filesize > 0) {
        pngdata <- .Call("b64encode", readBin(filename, "raw", n = filesize),
                         PACKAGE = "gridSVG")
        paste0("data:image/png;base64,", pngdata, collapse = "")
    } else {
        warning(paste(sQuote(filename), "is empty"))
        filename
    }
}
