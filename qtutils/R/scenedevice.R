##' Starts a QGraphicsScene graphics device.
##'
##' \code{qsceneDevice} starts a graphics device that is based on Qt's
##' Graphics View framework.  In the abstract sense, the device is a
##' \sQuote{graphics scene}, which contains various graphical items
##' such as circles, lines, and text.  The scene itself is not tied to
##' a particular on-screen view or output file.  Further steps must be
##' taken to view the contents of the scene, or render it to a file in
##' a suitable format.  See \code{\link{QT}} for a convenient wrapper
##' that provides these.
##' 
##' @title qsceneDevice
##' @param width Width of the scene in inches, assuming 72 pixels per inch.
##' @param height Height of the scene in inches, assuming 72 pixels per inch.
##' @param pointsize Default pointsize.
##' @param family Default font family.
##' 
##' @return A QGraphicsScene instance (same as the \code{rscene}
##' argument).  Drawing operations will result in QGraphicsItems being
##' added to the scene.  Note that unlike common R graphics devices,
##' the return value is nontrivial.
##' 
##' @author Deepayan Sarkar
qsceneDevice <-
    function(width = 10, height = 10, pointsize = 12, family = "")
    ## , rscene = Qt$QGraphicsScene()
    ##  @param rscene A QGraphicsScene instance.  If missing, a new
    ##                instance will be created.
{
    ## force(rscene)
    .Call(qt_qsceneDevice,
          as.numeric(width),
          as.numeric(height),
          as.numeric(pointsize),
          as.character(family))
    ## rscene)
    ## invisible(rscene)
}

##' Convenience wrapper for the \code{\link{qsceneDevice}} graphics device.
##'
##' \code{QT} is convenient user-level wrapper for the
##' \code{\link{qsceneDevice}} graphics device.  It returns a
##' QGraphicsView instance (which is also a QWidget instance) whose
##' scene is set to the QGraphicsScene instance created by a call to
##' \code{qsceneDevice(\dots)}.  In addition, several predefined
##' actions are associated with the view, allowing (through a context
##' menu and keyboard shortcuts) zooming, printing, and exporting as a
##' bitmap image.
##' 
##' @title QT
##' 
##' @param rscene A QGraphicsScene instance produced by a call to
##' \code{\link{qsceneDevice}}.  Can be missing, in which case a
##' suitable instance will be created (see \code{\dots} below).
##' @param ... Arguments passed on to \code{\link{qsceneDevice}} if
##' \code{rscene} is missing.
##' @param antialias Logical flag.  Specifies whether the view should
##' be antialiased.
##' @param opengl Logical flag.  Specifies whether the view should be
##' a QGLWidget, used for rendering OpenGL graphics.
##' 
##' @return A QGraphicsView instance
##' 
##' @author Deepayan Sarkar
QT <- function(rscene, ..., antialias = TRUE, opengl = FALSE)
{
    if (missing(rscene)) rscene <- qsceneDevice(...)
    gview <- Qt$QGraphicsView(rscene)
    if (antialias) gview$setRenderHints(Qt$QPainter$Antialiasing) 
    gview$setDragMode(Qt$QGraphicsView$ScrollHandDrag)

    ## Activate context menu with actions
    gview$setContextMenuPolicy(Qt$Qt$ActionsContextMenu)

    ## Add "Actions" to scale
    zoominAct <- Qt$QAction("Zoom In", gview)
    zoominAct$setShortcut(Qt$QKeySequence("Ctrl++"))
    qconnect(zoominAct,
             signal = "triggered",
             handler = function(checked) {
                 gview$scale(1.2, 1.2)
             })
    gview$addAction(zoominAct)
    zoomoutAct <- Qt$QAction("Zoom Out", gview)
    zoomoutAct$setShortcut(Qt$QKeySequence("Ctrl+-"))
    qconnect(zoomoutAct,
             signal = "triggered",
             handler = function(checked) {
                 gview$scale(1/1.2, 1/1.2)
             })
    gview$addAction(zoomoutAct)

    ## Helper function to print
    printHandler <- function(full = TRUE)
    {
        printer <- Qt$QPrinter(Qt$QPrinter$HighResolution)
        rpaper <- getOption("papersize")
        if (is.null(rpaper)) rpaper <- "A4"
        qtpaper <- names(Qt$QPrinter)
        usepaper <- qtpaper[ match(tolower(rpaper), tolower(qtpaper)) ]
        if (is.na(usepaper)) usepaper <- "A4"
        printer$setPageSize(Qt$QPrinter[[usepaper]])
        pd <- Qt$QPrintDialog(printer)
        acceptPrint <- pd$exec()
        if (acceptPrint)
        {
            painter <- Qt$QPainter()
            painter$begin(printer)
            if (full)
                rscene$render(painter)
            else
                gview$render(painter)
            painter$end()
        }
    }
    
    ## Actions to print
    printAct <- Qt$QAction("Print", gview)
    printAct$setShortcut(Qt$QKeySequence("Ctrl+P"))
    qconnect(printAct,
             signal = "triggered",
             handler = function(checked) {
                 printHandler(TRUE)
             })
    gview$addAction(printAct)
    printVisibleAct <- Qt$QAction("Print visible", gview)
    qconnect(printVisibleAct,
             signal = "triggered",
             handler = function(checked) {
                 printHandler(FALSE)
             })
    gview$addAction(printVisibleAct)

    ## Action to export (to image file).  This is a bit more involved,
    ## so factor out as separate function.
    addImageExportAction(gview)

    if (opengl)
      gview$setViewport(Qt$QGLWidget())

    gview$setWindowTitle("[ACTIVE] QGraphicsView (QGraphicsScene) Device")
    ## Return view widget
    gview
}

addImageExportAction <- function(gview)
{
    saveAsMenu <- Qt$QMenu("&Export As", gview)
    supportedFormats <-
        unique(toupper(sapply(Qt$QImageWriter$supportedImageFormats(),
                              rawToChar)))
    for (fmt in supportedFormats)
    {
        action <- Qt$QAction(fmt, gview)
        qconnect(action, signal = "triggered",
                 handler = function(checked, user.data) {
                     saveAsImage(gview, user.data)
                 }, user.data = force(fmt))
        saveAsMenu$addAction(action)
    }
    saveAsAct <- saveAsMenu$menuAction()
    gview$addAction(saveAsAct)
}

saveAsImage <- function(gview, format, full = TRUE)
{
    initialPath <- file.path(getwd(), paste("untitled.", tolower(format), sep = ""))
    fileName <-
        Qt$QFileDialog$getSaveFileName(gview, "Save As",
                                       initialPath,
                                       sprintf("%s Files (*.%s);;All Files (*)",
                                               toupper(format),
                                               toupper(format)))
    if (!is.null(fileName)) {
        qexport(gview, file = fileName, format = format, full = full)
    }
}



