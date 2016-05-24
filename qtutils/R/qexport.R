
##' Exports a Qt object, usually as a bitmap image.
##'
##' A generic function that is meant to export a Qt object to an
##' external file.  Currently, only export to some image formats is
##' supported.  This is used, in particular, to allow exporting the
##' contents of a \code{\link{qsceneDevice}} via a context menu item
##' for a corresponding view. 
##' 
##' @title Export Qt object
##' @param x 
##' @param ... Passed on to methods
##' @author Deepayan Sarkar
qexport <- function(x, ...)
{
    UseMethod("qexport")
}

qexport.QGraphicsScene <-
    function(x, file,
             format = tail(strsplit(basename(file), ".", fixed = TRUE)[[1]], 1),
             bg, ...)
{
    if (missing(bg))
        bg <- switch(toupper(format),
                     JPG=, JPEG=, BMP= "white",
                     "transparent")
    size <- x$sceneRect$size()
    qimg <- Qt$QImage(size$toSize(), Qt$QImage$Format_ARGB32_Premultiplied)
    painter <- Qt$QPainter()
    painter$begin(qimg)
    if (bg != "transparent") painter$fillRect(x$sceneRect, qcolor(bg))
    painter$setRenderHint(Qt$QPainter$Antialiasing)
    painter$setRenderHint(Qt$QPainter$TextAntialiasing)
    x$render(painter)
    painter$end()
    qimg$save(file, toupper(format))
}

qexport.QWidget <- 
    function(x, file,
             format = tail(strsplit(basename(file), ".", fixed = TRUE)[[1]], 1),
             bg = "white", ...)
{
    size <- x$size
    qimg <- Qt$QImage(size$toSize(), Qt$QImage$Format_ARGB32_Premultiplied)
    painter <- Qt$QPainter()
    painter$begin(qimg)
    if (bg != "transparent") painter$fillRect(x$rect, qcolor(bg))
    painter$setRenderHint(Qt$QPainter$Antialiasing)
    painter$setRenderHint(Qt$QPainter$TextAntialiasing)
    x$render(painter)
    painter$end()
    qimg$save(file, toupper(format))
}


qexport.QGraphicsView <- function(x, ..., full = FALSE)
{
    ## Is this doing what we mean to do?
    if (full) qexport(x$scene(), ...)
    else NextMethod("qexport")
}



## May want to do something common for printing, and maybe also
## PDF/SVG export to file with custom width/height and margins.
## Direct export to PDF (like copy2pdf) could definitely be useful.


## qexport.QWidget <-
##     function(x,
##              file = qfile.choose(caption = "Choose output file",
##                                  allow.new = TRUE, parent = x),
##              type, ...)
## {
##     file <- as.character(file)
##     stopifnot(length(file) == 1)
##     if (nzchar(file)) ## file=="" if selection cancelled
##     {
##         extension <- tail(strsplit(basename(file), ".", fixed = TRUE)[[1]], 1)
##         if (missing(type))
##             type <-
##                 if (tolower(extension) %in% c("ps", "pdf")) "vector"
##                 else if (tolower(extension) %in% "svg") "svg"
##                 else "raster"
##         switch(type,
##                svg = qrenderToSVG(x, file),
##                raster = qrenderToPixmap(x, file),
##                vector = qrender(x, file))
##         invisible()
##     }
## }
