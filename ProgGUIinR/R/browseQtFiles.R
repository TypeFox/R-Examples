##' @include misc.R
NULL

##' Creates GUI to browse QT examples
##' 
##' @return makes the GUI to browse the example files for Qt.
##' @export
browseQtFiles <- function() {
  if(!faux_require("qtbase"))
    stop("This function needs the qtbase package")

  browseFiles <- function() {}
  source(system.file("qt", "browseQtFiles.R", package="ProgGUIinR"), local=FALSE)
  
  w <- browseFiles()
  w$show()
  invisible(w)
}


          
