##' @include misc.R
NULL

##' Simple GUI to browse gWidgets examples
##' 
##' @param toolkit Which toolkit to use (RGtk2, tcltk or Qt)
##' @return NULL
##' @export
browsegWidgetsFiles <- function(toolkit="RGtk2") {
  e <- new.env()                        # for evaluation

  if(!faux_require("gWidgets"))
    stop("Requires gWidgets package")

  options(guiToolkit=toolkit)
  
  
  listFiles <- function() {
    files <- list.files(path=dir, pattern="^ex-.*R$", full.names=TRUE)
    names(files) <- list.files(path=dir, pattern="^ex-.*R$", full.names=FALSE)
    return(files)
  }

  dir <- system.file("Examples","ch-gWidgets", package="ProgGUIinR")
  files <- listFiles()
  
  w <- gwindow("Browse files", width=800, height=500, visible=FALSE)
  gstatusbar("Select a file to browse or run its source code", container=w)
  
  pg <- gpanedgroup(container = w, expand=TRUE)
  f <- gframe("gWidgets Examples", container=pg)
  fBrowser <- gtable(data.frame(Files=names(files), stringsAsFactors=FALSE), expand=TRUE, container = f,
                     handler = function(h, ...) {
                       val <- svalue(h$obj)
                       f <- files[val]
                       nm <- gsub("ex-gWidgets-","",val)
                       if(nm %in% names(sourceNb)) {
                         svalue(sourceNb) <- which(nm == names(sourceNb))
                       } else {
                         t <- gtext("", container=sourceNb, label=nm)
                         svalue(t) <- readLines(f)
                         svalue(sourceNb) <- length(sourceNb) # set tab
                       }
                     })

  g <- ggroup(horizontal=FALSE, container = pg)

  bg <- ggroup(container = g, horizontal=TRUE); addSpring(bg)
  runButton <- gbutton("source", container = bg, handler = function(h,...) {
    if(length(sourceNb) == 0) return()
    t <- sourceNb[svalue(sourceNb)]     # get child
    out <- paste(svalue(t), collapse="\n")
    eval(parse(text=out), envir=e)
  })

  sourceNb<- gnotebook(container=g, expand=TRUE, closebuttons=TRUE)

  visible(w) <- TRUE
  svalue(pg) <- 0.33                    # adjust sash
}

