##' @title Start the Retistruct GUI
##' @seealso gWidgets
##' @return Object with \code{getData()} method to return
##' reconstructed retina data and environment \code{this} which
##' contains variables in object.
##' @export
retistruct <- function() {
  ## This function is essentially the constructor for a class. The
  ## environment 'this' contains all member variables and functions of
  ## the class. We do not need to refer to 'this' in the code, but we
  ## will return it to facilitate debugging. 
  this <- environment()

  ## Global variables
  dataset <- NULL                         # Directory of dataset
  initial.dir <- "."
  a <- NULL                             # Annotation object
  r <- NULL                             # Reconstruction object
  ## Path to extdata for demos
  extdata       <- file.path(system.file(package = "retistruct"), "extdata")
  extdata.demos <- file.path(system.file(package = "retistructdemos"), "extdata")
  
  ## Accessor functions
  getData <- function() {
    return(r)
  }

  ##
  ## Load and install graphics packages
  ##
  
  ## TODO: It would be good to make it possible to set
  ## guiToolkit=tcltk as an option to the function so that there is no
  ## need to install gWidgetsRGtk2, which is problematic on a Mac
  ##
  ## @param guiToolkit The toolkit \code{gWidgets} toolkit that will
  ## be ' used for the GUI
  guiToolkit <- "RGtk2"
  require.package <- function(pkg) {
    if (!require(pkg, character.only=TRUE)) {
      message(paste("Trying to install required package", pkg))
      install.packages(pkg)
      if (!require(pkg, character.only=TRUE)) {
        stop(paste("Could not install", pkg))
      }
    }
  }
  require.package(paste0("gWidgets", guiToolkit))
  require.package("cairoDevice")
  options(guiToolkit=guiToolkit)
  
  ##
  ## Convenience functions for handlers
  ## 

  enable.group <- function(widgets, state=TRUE) {
    for (w in widgets) {
      gWidgets::enabled(w) <- state
    }
  }

  enable.widgets <- function(state) {
    enable.group(c(g.add, g.move, g.remove, g.reconstruct, g.properties,
                   g.mark.n, g.mark.d, g.mark.od,
                   g.phi0d, g.show, g.edit.show, g.data, g.eye,
                   g.print1, g.print2,
                   g.print.pdf1, g.print.pdf2), state)
    if (state) 
      enable.group(c(g.mark.od), retistruct.potential.od(a))
    if (!retistruct.check.markup(a)) {
      enable.group(c(g.reconstruct), FALSE)
    }
  }

  unsaved.data <- function(state) {
    if (state) {
      r <<- NULL
    }
    enable.group(c(g.save), state)
  }

  ## Function to report to set status
  set.status <- function(s) {
    gWidgets::svalue(g.status) <- s
  }

  ## Utility function
  identify.abort.text <- function() {
    if (.Platform$GUI == "X11") {
      return("Right-click to abort.")
    }
    if (.Platform$GUI == "AQUA") {
      return("Press ESC to abort.")
    }
  }

  ##
  ## Editting handlers
  ##

  ## Handler for adding a point
  h.add <- function(h, ...) {
    unsaved.data(TRUE)
    enable.widgets(FALSE)
    gWidgets::svalue(g.status) <- paste("Click on the three points of the tear in any order.",
                              identify.abort.text())
    dev.set(d1)
    pids <- with(a, identify(P[,1], P[,2], n=3, col=getOption("TF.col")))
    withCallingHandlers({
      a <<- addTear(a, pids)
    }, warning=h.warning, error=h.warning)
    do.plot()
    gWidgets::svalue(g.status) <- ""
    enable.widgets(TRUE)
  }

  ## Handler for removing points
  h.remove <- function(h, ...) {
    unsaved.data(TRUE)
    enable.widgets(FALSE)
    gWidgets::svalue(g.status) <- paste("Click on the apex of the tear to remvoe.",
                              identify.abort.text())
    dev.set(d1)
    id <- with(a, identify(P[,1], P[,2], n=1, plot=FALSE))
    a <<- removeTear(a, whichTear(a, id))
    do.plot()
    gWidgets::svalue(g.status) <- ""
    enable.widgets(TRUE)
  }

  ## Handler for moving points
  h.move <- function(h, ...) {
    unsaved.data(TRUE)
    enable.widgets(FALSE)
    dev.set(d1)
    ## Find the intial point
    gWidgets::svalue(g.status) <- paste("Click on apex or vertex to move.",
                              identify.abort.text())
    id1 <- with(a, identify(P[,1], P[,2], n=1, plot=FALSE))
    
    ## Locate tear ID in which the point occurs
    tid <- whichTear(a, id1)

    ## If there is a tear in which it occurs, select a point to move it to
    if (!is.na(tid)) {
      gWidgets::svalue(g.status) <- paste("Click on point to move it to.",
                                identify.abort.text())

      ## Label first point
      with(a, points(P[id1,1], P[id1,2], col="yellow"))

      ## Select second point
      id2 <- with(a, identify(P[,1], P[,2], n=1))

      ## Get point ids of exsiting tear
      pids <- getTear(a, tid)

      ## Replace old point with desired new point
      if (length(id2)) pids[pids==id1] <- id2

      ## It is possible to get the apex and vertex mixed up when moving points.
      ## Fix any errors.
      pids <- labelTearPoints(a, pids)
      a <<- removeTear(a, tid)
      a <<- addTear(a, pids)
    }

    ## Display and cleanup
    do.plot()
    gWidgets::svalue(g.status) <- ""
    enable.widgets(TRUE)
  }

  ## Handler for marking nasal point
  h.mark.n <- function(h, ...) {
    unsaved.data(TRUE)
    enable.widgets(FALSE)
    gWidgets::svalue(g.status) <- paste("Click on nasal point.",
                              identify.abort.text())
    dev.set(d1)
    id <- with(a, identify(P[,1], P[,2], n=1))
    withCallingHandlers({
      a <<- setFixedPoint(a, id, "Nasal")
    }, warning=h.warning, error=h.warning)
    do.plot()
    gWidgets::svalue(g.status) <- ""
    enable.widgets(TRUE)
  }

  ## Handler for marking dorsal point
  h.mark.d <- function(h, ...) {
    unsaved.data(TRUE)
    enable.widgets(FALSE)
    gWidgets::svalue(g.status) <- paste("Click on dorsal point.",
                              identify.abort.text())
    dev.set(d1)
    id <- with(a, identify(P[,1], P[,2], n=1))
    withCallingHandlers({
      a <<- setFixedPoint(a, id, "Dorsal")
    }, warning=h.warning, error=h.warning)
    do.plot()
    gWidgets::svalue(g.status) <- ""
    enable.widgets(TRUE)
  }

  ## Handler for marking optic disc
  h.mark.od <- function(h, ...) {
    unsaved.data(TRUE)
    enable.widgets(FALSE)
    gWidgets::svalue(g.status) <- paste("Click on a point on the optic disc.",
                              identify.abort.text())
    dev.set(d1)
    ## Convert list of segments to a matrix of points
    Sm <- NULL
    for (S in a$Ss) {
      Sm <- rbind(Sm, S)
    }

    ## Identify a point
    id <- identify(Sm[,1], Sm[,2], n=1)
    
    ## Idendify segment in which point appears
    i <- 0
    N <- 0
    while (id > N && i < length(a$Ss)) {
      i <- i + 1
      N <- N +  nrow(a$Ss[[i]])
    }
    ## Set "OD" landmark
    a <<- nameLandmark(a, i, "OD")
    do.plot()
    gWidgets::svalue(g.status) <- ""
    enable.widgets(TRUE)
  }

  ## Handler for setting phi0d
  h.phi0d <- function(h, ...) {
    unsaved.data(TRUE)
    v <- gWidgets::svalue(g.phi0d)
    if (v < -80) {
      v <- -89
    }
    if (v > 89) {
      v <- 89
    }
    a$phi0 <<- v*pi/180
  }

  ## Handler for saving state
  h.save <- function(h, ...) {
    retistruct.save.markup(a)
    ## If the reconstruction doesn't exist, remove the reconstruction
    ## file to ensure consistency
    if (is.null(r)) {
      unlink(file.path(a$dataset, "r.Rdata"))
    } else {
      retistruct.save.recdata(r)
    }
    retistruct.export.matlab(r)
    unsaved.data(FALSE)
  }

  ## Handler for brining up a file dialogue to open a dataset
  ##
  ## Changes the global object r
  ##
  ## Produces a plot of the retina in device d1
  ## 
  h.select <- function(h, ...) {
    curdir <- getwd()
    if (is.null(a$dataset)) {
      info = file.info(initial.dir)
      if (!is.na(info$isdir)) {
        setwd(initial.dir)
      }
    } else {
      setwd(a$dataset)
      setwd("..")
    } 
    gWidgets::gfile(type="selectdir", text="Select a directory...",
          handler = function(h, ...) {
            a$dataset <<- h$file
          })
    setwd(curdir)
    h.open()
  }

  ## Handler for opening a file
  h.open <- function(h, ...) {
    ## Read the raw data
    withCallingHandlers({
      a <<- retistruct.read.dataset(a$dataset)
    }, warning=h.warning, error=h.error)
    
    ## Read the markup
    withCallingHandlers({
      a <<- retistruct.read.markup(a, error=message)
    }, warning=h.warning, error=h.warning)
    gWidgets::svalue(g.win)   <- paste(version.string(), "-" ,a$dataset)
    gWidgets::svalue(g.phi0d) <- a$phi0*180/pi
    gWidgets::svalue(g.eye)   <- a$side
    
    ## Read the reconstruction data
    withCallingHandlers({
      r <<- retistruct.read.recdata(a, check=TRUE)
    }, warning=h.warning, error=h.error)
    ## If there is no reconstruction data, show the markup so that we
    ## don't think there is no markup.
    if (is.null(r)) {
      gWidgets::svalue(g.show) <- unique(c("Markup", gWidgets::svalue(g.show)))
      gWidgets::svalue(g.nb) <- 1                   # Set "Edit" tab
    } else {
      gWidgets::svalue(g.nb) <- 2                   # Set "View" tab
    }
    gWidgets::delete(g.ids.frame, g.ids)
    ids <- getIDs(a)
    if (!is.null(ids)) {
      g.ids <<- gWidgets::gcheckboxgroup(ids, checked=rep(TRUE, length(ids)),
                               handler=h.show, container=g.ids.frame)
    }
    unsaved.data(FALSE)
    enable.widgets(TRUE)
    do.plot()
  }

  ## Handler to start reconstructing the retina
  h.reconstruct <- function(h, ...) {
    unsaved.data(TRUE)
    enable.widgets(FALSE)
    withCallingHandlers({
      r <<- retistruct.reconstruct(a, report=set.status,
                                   plot.3d=getOption("show.sphere"),
                                   dev.flat=d1, dev.polar=d2)
    }, warning=h.warning, error=h.warning)  
    enable.widgets(TRUE)
    do.plot()
  }

  ## Handler for showing data
  h.show <- function(h, ...) {
    if(!is.null(h$pageno)) {
      if (h$pageno == 1) {
        do.plot(markup=TRUE)
      } else {
        do.plot(markup=("Markup" %in% (gWidgets::svalue(g.show))))
      }
    } else {
      do.plot()
    }
  }

  ## Handler for flipping DV axis
  h.flipdv <- function(h, ...) {
    unsaved.data(TRUE)
    a$DVflip <<- ("Flip DV" %in% gWidgets::svalue(g.data))
    do.plot()
  }

  ## Handler for dealing with data
  h.eye <- function(h, ...) {
    unsaved.data(TRUE)
    a$side <<- gWidgets::svalue(g.eye)
    do.plot()
  }

  ## Print device d to file
  print.bitmap <- function(d, file) {
    dev <- NULL
    if (grepl(".png$", file, ignore.case=TRUE)) 
      dev <- png
    if (grepl(".jpeg$", file, ignore.case=TRUE) ||
        grepl(".jpg$", file, ignore.case=TRUE))
      dev <- jpeg
    if (grepl(".tif$", file, ignore.case=TRUE) ||
        grepl(".tiff$", file, ignore.case=TRUE))
      dev <- tiff
    if (is.null(dev)) {
      file <- paste(file, ".png")
      dev <- png
    }
    dev.set(d)
    dev.print(dev, file, width=1000, height=1000)
  }

  ## Print device d to file
  print.pdf <- function(d, file) {
    dev.set(d)
    dev.print(pdf, file,
              width=getOption("retistruct.print.pdf.width"),
              height=getOption("retistruct.print.pdf.width"))
  }

  ## Handler for printing to a bitmap
  h.print.bitmap <- function(d, initialfilename) {
    curdir <- getwd()
    setwd(a$dataset)  
    gWidgets::gfile(type="save", text="Select a filename to save image to...",
          initialfilename=initialfilename,
          handler=function(h, ...) {
            print.bitmap(d, h$file)
          })
    setwd(curdir)
  }

  ## Handler for printing to a pdf file
  h.print.pdf <- function(d, initialfilename) {
    h.width <- function(h, ...) {
      options(retistruct.print.pdf.width=gWidgets::svalue(g.width))
    }
    g.pdf <- gWidgets::gbasicdialog(title="PDF options", 
                          handler=h.width)
    g.pdf.group <- gWidgets::ggroup(container=g.pdf, horizontal=TRUE)
    gWidgets::glabel("Width & height (inches)", container=g.pdf.group)
    g.width <- gWidgets::gedit(getOption("retistruct.print.pdf.width"),
                     width=5, coerce.with=as.numeric,
                     container=g.pdf.group)


    gWidgets::visible(g.pdf, set=TRUE)

    curdir <- getwd()
    setwd(a$dataset)  
    gWidgets::gfile(type="save", text="Select a filename to save image to...",
          initialfilename=initialfilename,
          handler=function(h, ...) {
            print.pdf(d, h$file)
          })
    setwd(curdir)
  }

  ## Handlers for printing bitmaps and PDFs from the two graphics
  ## devices
  h.print1 <- function(h, ...) {
    h.print.bitmap(d1, initialfilename="image-flat.png")
  }

  h.print.pdf1 <- function(h, ...) {
    h.print.pdf(d1, initialfilename="image-flat.pdf")
  }

  h.print2 <- function(h, ...) {
    h.print.bitmap(d2, initialfilename="image-polar.png")
  }

  h.print.pdf2 <- function(h, ...) {
    h.print.pdf(d2, initialfilename="image-polar.pdf")
  }

  ## Get and name available projections
  getProjections <- function() {
    return(list("Azimuthal Equidistant"=azimuthal.equidistant,
                "Azimuthal Equal Area" =azimuthal.equalarea,
                "Azimuthal Conformal"  =azimuthal.conformal,
                "Sinusoidal"           =sinusoidal,
                "Orthographic"         =orthographic))
  }

  ## Get and name available transforms
  getTransforms <- function() {
    return(list("None"  =identity.transform,
                "Invert"=invert.sphere,
                "Invert to hemisphere"=invert.sphere.to.hemisphere))
  }

  ## Plot in edit pane
  do.plot <- function(markup=("Markup" %in% (gWidgets::svalue(g.show))) | (gWidgets::svalue(g.nb) == 1)) {
    
    if (is.null(r)) {
      r <- a
    }
    if ("Strain" %in% gWidgets::svalue(g.edit.show)) {   # Strain plot
      dev.set(d1)
      par(mar=c(0.5, 0.5, 0.5, 0.5))
      flatplot(r, axt="n",
               datapoints=FALSE,
               landmarks=FALSE,
               markup=FALSE,
               stitch=FALSE,
               grid=FALSE,
               mesh=FALSE,
               strain=TRUE,
               scalebar=1)
      dev.set(d2)
      par(mar=c(4.5, 4.5, 1, 0.5))
      lvsLplot(r)
    } else {
      dev.set(d1)
      par(mar=c(0.5, 0.5, 0.5, 0.5))
      flatplot(r, axt="n",
               datapoints=("Points" %in% gWidgets::svalue(g.show)),
               grouped=("Counts" %in% gWidgets::svalue(g.show)),
               landmarks=("Landmarks" %in% gWidgets::svalue(g.show)),
               markup=markup,
               stitch=("Stitch" %in% gWidgets::svalue(g.show)),
               grid=("Grid" %in% gWidgets::svalue(g.show)),
               ids=gWidgets::svalue(g.ids),
               mesh=FALSE,
               scalebar=1)
      dev.set(d2)
      par(mar=c(0.7, 0.7, 0.7, 0.7))
      projection(r,
                 datapoints=("Points" %in% gWidgets::svalue(g.show)),
                 datapoint.means=("Point means" %in% gWidgets::svalue(g.show)),
                 landmarks=("Landmarks" %in% gWidgets::svalue(g.show)),
                 transform=getTransforms()[[gWidgets::svalue(g.transform)]],
                 projection=getProjections()[[gWidgets::svalue(g.projection)]],
                 axisdir=cbind(phi=gWidgets::svalue(g.axis.el), lambda=gWidgets::svalue(g.axis.az)),
                 proj.centre=cbind(phi=gWidgets::svalue(g.pc.el), lambda=gWidgets::svalue(g.pc.az)),
                 datapoint.contours=("Point contours" %in% gWidgets::svalue(g.show)),
                 grouped=("Counts" %in% gWidgets::svalue(g.show)),
                 grouped.contours=("Count contours" %in% gWidgets::svalue(g.show)),
                 ids=gWidgets::svalue(g.ids))
      ## FIXME: EOD not computed
      if (!is.null(r$EOD)) {
        polartext(paste("OD displacement:",
                        format(r$EOD, digits=3, nsmall=2), "deg"))
      }
    }
    dev.set(d1)
  }

  ## It would be nice to have error messages displayed graphically.
  ## This function should work, but has the problem that it always gives an
  ## "Error in get(\"toolkit\", inherits = TRUE) : object 'toolkit' not found\n"
  ## error itself.
  h.error <- function(e) {
    gWidgets::gmessage(e, title="Error", icon="error")
    stop(e)
  }

  ## Warning message
  h.warning <- function(e) {
    gWidgets::gmessage(e, title="Warning", icon="warning")
    invokeRestart("muffleWarning")
  }


  ## Poperties dialogue
  h.properties <- function(h, ...) {
    g.win <- gWidgets::gwindow("Properties",
                     parent=g.win)
    g.props <- gWidgets::ggroup(container=g.win, horizontal=FALSE)
    g.colours <- gWidgets::gframe("Colours", container=g.props, horizontal=FALSE)
    cols <- c("black", "red", "green3", "blue", "cyan", "magenta", "yellow",
              "gray")
    g.prop.dl <- function(name, property, container) {
      g.prop.dl.group <- gWidgets::ggroup(container=container)
      gWidgets::glabel(name, container=g.prop.dl.group)
      g.dl <- gWidgets::gdroplist(cols,
                        selected=which(cols == options(property[1])),
                        container=g.prop.dl.group,
                        handler=function(h, ...) {
                          for (p in property) {
                            eval(parse(text=paste("options(", p, "=gWidgets::svalue(g.dl))")))
                          }
                          do.plot()})
    }

    g.prop.dl("Outline colour", "outline.col", g.colours)
    g.prop.dl("Tear colour", c("TF.col", "TB.col", "V.col"), g.colours)
    g.prop.dl("Stitch colour", "stitch.col", g.colours)
    g.prop.dl("Major gridline colour", "grid.maj.col", g.colours)
    g.prop.dl("Minor gridline colour", "grid.min.col", g.colours)

    g.printing <- gWidgets::gframe("Printing", container=g.props, horizontal=FALSE)
    g.group <- gWidgets::ggroup(container=g.printing)
    gWidgets::glabel("Maximum resolution of projection", container=g.group)
    property <- "max.proj.dim"
    g.max.proj.dim <- gWidgets::gedit(0, width=5, coerce.with=as.numeric,
                            container=g.group,
                            handler=function(h, ...) {
                              eval(parse(text=paste("options(", property, "=gWidgets::svalue(g.max.proj.dim))")))
                            })
    gWidgets::svalue(g.max.proj.dim) <- getOption(property)
    
    gWidgets::gbutton("Close", container=g.props,
            handler = function(h,...) gWidgets::dispose(g.win))
  }

  ## Construct the version string
  version.string <- function() {
    return(paste0("Retistruct ",
                  packageDescription("retistruct", fields="Version"),
                  " (Revision", retistruct.global.revision, " of ",
                  packageDescription("retistruct", fields="Date"), ")"))
  }
  
  ##
  ## GUI Layout
  ##

  ## Top level window
  g.win <- gWidgets::gwindow(version.string())

  ## Menu in row 0
  mbl <- list()
  mbl$File$Open$handler <- h.select
  mbl$File$Save$handler <- h.save
  mbl$Edit$Reconstruct$handler <- h.reconstruct
  mbl$Edit$Properties$handler <- h.properties
  mbl$Edit$Properties$handler <- h.properties
  mbl$Demos$fig1 <-
    gWidgets::gaction(label="Figure 1",
            handler=function(h, ...) {
              a$dataset <<- file.path(extdata, "GM509", "R-CONTRA")
              h.open()
            })
  mbl$Demos$fig2low <-
    gWidgets::gaction(label="Figure 2A-D: Low deformation",
            handler=function(h, ...) {
              a$dataset <<- file.path(extdata, "GMB530", "R-CONTRA")
              h.open()
            })
  mbl$Demos$fig2high <-
    gWidgets::gaction(label="Figure 2E-H: High deformation",
            handler=function(h, ...) {
              a$dataset <<- file.path(extdata, "GM182-4", "R-CONTRA")
              h.open()
            })
  mbl$Demos$SMI32$handler <- function(h, ...) {
    a$dataset <<- file.path(extdata, "smi32")
    h.open()
  }
  mbl$Demos$left.contra <-
    gWidgets::gaction(label="Figure 6 Left Contra",
            handler=function(h, ...) {
              a$dataset <<- file.path(extdata.demos, "Figure_6-data", "left-contra")
              h.open()
            })
  mbl$Demos$left.ipsi <-
    gWidgets::gaction(label="Figure 6 Left Ipsi",
            handler=function(h, ...) {
              a$dataset <<- file.path(extdata.demos, "Figure_6-data", "left-ipsi")
              h.open()
            })
  mbl$Demos$right.contra <-
    gWidgets::gaction(label="Figure 6 Right Contra",
            handler=function(h, ...) {
              a$dataset <<- file.path(extdata.demos, "Figure_6-data", "right-contra")
              h.open()
            })
  mbl$Demos$right.ipsi <-
    gWidgets::gaction(label="Figure 6 Right Ipsi",
            handler=function(h, ...) {
              a$dataset <<- file.path(extdata.demos, "Figure_6-data", "right-ipsi")
              h.open()
            })
  mbl$Help$About$handler <- function(h, ...) {
    gWidgets::gmessage("Retistruct was written by David Sterratt at the University of Edinburgh, and tested by Daniel Lyngholm and Ian Thompson at the MRC Centre for Developmental Neurobiology, KCL.

This work was supported by a Programme Grant from the Wellcome Trust (G083305). ", title="About")
  }

  g.menu <- gWidgets::gmenu(mbl, container=g.win)
  
  g.rows <- gWidgets::ggroup(horizontal=FALSE, container=g.win)
  ## Toolbar in row 1
  g.open         <- gWidgets::gaction("Open", icon="open", handler=h.select)
  g.save         <- gWidgets::gaction("Save", icon="save", handler=h.save)
  g.reconstruct  <- gWidgets::gaction("Reconstuct retina", icon="polar", handler=h.reconstruct)
  g.properties   <- gWidgets::gaction("Properties", icon="properties", handler=h.properties)
  g.toolbar <- gWidgets::gtoolbar(list(open=g.open,
                             save=g.save,
                             reconstruct=g.reconstruct,
                             options=g.properties),
                        container=g.rows, style="both")

  ## Body of interface
  g.body <- gWidgets::ggroup(container=g.rows)

  ## "Edit" and "View" tabs
  g.nb <- gWidgets::gnotebook(container=g.body)

  ## Edit tab
  
  ## Tear editor
  g.editor <- gWidgets::ggroup(horizontal = FALSE, container=g.nb, label="Edit")

  g.add     <- gWidgets::gbutton("Add tear",    handler=h.add,     container=g.editor)
  g.move    <- gWidgets::gbutton("Move Point",  handler=h.move,    container=g.editor)
  g.remove  <- gWidgets::gbutton("Remove tear", handler=h.remove,  container=g.editor)
  g.mark.n  <- gWidgets::gbutton("Mark nasal",  handler=h.mark.n,  container=g.editor)
  g.mark.d  <- gWidgets::gbutton("Mark dorsal", handler=h.mark.d,  container=g.editor)
  g.mark.od <- gWidgets::gbutton("Mark OD",     handler=h.mark.od, container=g.editor)
  
  ## Editting of data
  g.data.frame <- gWidgets::gframe("Data", container=g.editor, horizontal=FALSE)
  g.data <- gWidgets::gcheckboxgroup(c("Flip DV"),
                           checked=c(FALSE),
                           handler=h.flipdv, container=g.data.frame)
  g.eye.frame <- gWidgets::gframe("Eye", container=g.editor, horizontal=FALSE)
  g.eye <- gWidgets::gradio(c("Right", "Left"),
                  checked=c(FALSE),
                  handler=h.eye, container=g.eye.frame)
  
  ## Editing of phi0
  g.phi0d.frame <- gWidgets::gframe("Phi0", container=g.editor)
  g.phi0d <- gWidgets::gedit(0, handler=h.phi0d, width=5, coerce.with=as.numeric,
                   container=g.phi0d.frame)

  ## Whether to show strain
  g.edit.show.frame <- gWidgets::gframe("Show", container=g.editor)
  g.edit.show <- gWidgets::gcheckboxgroup(c("Strain"),
                                checked=c(FALSE),
                                handler=h.show, container=g.edit.show.frame)

  ## View Tab

  ## What to show
  g.view <- gWidgets::ggroup(horizontal=FALSE, container=g.nb, label="View")
  g.show.frame <- gWidgets::gframe("Show", container=g.view)
  g.show <- gWidgets::gcheckboxgroup(c("Markup", "Stitch", "Grid", "Landmarks",
                             "Points", "Point means", "Point contours",
                             "Counts", "Count contours"),
                           checked=c(TRUE, FALSE, FALSE, FALSE,
                             FALSE, FALSE, FALSE,
                             FALSE, FALSE),
                           handler=h.show, container=g.show.frame)

  ## Group IDs
  g.ids.frame <- gWidgets::gframe("IDs", container=g.view)
  g.ids <- gWidgets::gcheckboxgroup("All", checked=TRUE,
                          handler=h.show, container=g.ids.frame)

  
  ## Projection type
  g.projection.frame <- gWidgets::gframe("Projection", container=g.view)
  g.projection <- gWidgets::gdroplist(names(getProjections()), selected=1,
                            handler=h.show, 
                            action=NULL, container=g.projection.frame)

  ## Projection centre
  g.pc.frame <- gWidgets::gframe("Projection centre", container=g.view, horizontal=TRUE)
  gWidgets::glabel("El", container=g.pc.frame)
  g.pc.el <- gWidgets::gedit("0", handler=h.show, width=5, coerce.with=as.numeric,
                   container=g.pc.frame)
  gWidgets::glabel("Az", container=g.pc.frame)
  g.pc.az <- gWidgets::gedit("0", handler=h.show, width=5, coerce.with=as.numeric,
                   container=g.pc.frame)

  ## Transform
  g.transform.frame <- gWidgets::gframe("Transform", container=g.view)
  g.transform <- gWidgets::gdroplist(names(getTransforms()), selected = 1,  handler = h.show, 
                           action = NULL, container = g.transform.frame)

  ## Axis direction
  g.axisdir.frame <- gWidgets::gframe("Axis direction", container=g.view, horizontal=TRUE)
  gWidgets::glabel("El", container=g.axisdir.frame)
  g.axis.el <- gWidgets::gedit("90", handler=h.show, width=5, coerce.with=as.numeric,
                     container=g.axisdir.frame)
  gWidgets::glabel("Az", container=g.axisdir.frame)
  g.axis.az <- gWidgets::gedit("0", handler=h.show, width=5, coerce.with=as.numeric,
                     container=g.axisdir.frame)

  ## Graphs at right

  ## Flat plot
  g.f1 <- gWidgets::ggroup(horizontal=FALSE, container=g.body)
  ## Buttons
  g.f1.buttons <- gWidgets::ggroup(horizontal=TRUE, container=g.f1)
  g.print1     <- gWidgets::gbutton("Bitmap", handler=h.print1,     container=g.f1.buttons)
  g.print.pdf1 <- gWidgets::gbutton("PDF",    handler=h.print.pdf1, container=g.f1.buttons)
  ## Device itself
  g.fd1 <- gWidgets::ggraphics(expand=TRUE, width=500, height=500, ps=11, container=g.f1)
  d1 <- dev.cur()

  ## Projection
  g.f2 <- gWidgets::ggroup(horizontal=FALSE, container=g.body)
  ## Buttons  
  g.f2.buttons <- gWidgets::ggroup(horizontal=TRUE, container=g.f2)  
  g.print2     <- gWidgets::gbutton("Bitmap", handler=h.print2,     container=g.f2.buttons)
  g.print.pdf2 <- gWidgets::gbutton("PDF",    handler=h.print.pdf2, container=g.f2.buttons)
  ## Device itself
  g.fd2 <- gWidgets::ggraphics(expand=TRUE, width=500, height=500, ps=11, container=g.f2)
  d2 <- dev.cur()
  
  ## Status bar
  ## g.statusbar <- ggroup(container=g.rows)
  g.statusbar <- gWidgets::gframe("", expand=TRUE, container=g.rows)
  g.status <- gWidgets::glabel("", container=g.statusbar)
  gWidgets::addSpring(g.statusbar)

  ## Disable buttons initially
  unsaved.data(FALSE)
  enable.widgets(FALSE)

  ## Have to add the hander to the notebook at the end, otherwise
  ## there are complaints about various components not being defined.
  gWidgets::addHandlerChanged(g.nb, handler=h.show)

  return(invisible(list(getData=getData, env=this)))
}

