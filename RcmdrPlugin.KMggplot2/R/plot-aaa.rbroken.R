#' Base Class for Plot Windows
#'
#' \code{plot_base} class is the base class for plot windows.
#'
#' This class is a framework for implementing subclasses which show dialog boxes for graphics editing.
#'
#' @section Fields:
#' \describe{
#' \item{\code{top}: }{\code{tkwin} class object; parent of widget window.} 
#' \item{\code{alternateFrame}: }{\code{tkwin} class object; a special frame for some GUI parts.} 
#' \item{\code{rmlist}: }{List of character; deletable temporary objects.} 
#' \item{\code{mode}: }{numeric; the executive mode (0 = justDoIt, 1 = doItAndPrint).} 
#' }
#' @section Contains:
#' NULL
#' @section Methods:
#' \describe{
#' \item{\code{plotWindow()}: }{Create the window that make plots.}
#' \item{\code{savePlot(plot)}: }{Save the plot.}
#' \item{\code{registRmlist(object)}: }{Register deletable temporary objects.}
#' \item{\code{removeRmlist()}: }{Remove registered temporary objects.}
#' \item{\code{setFront()}: }{Set front parts of frames.}
#' \item{\code{setBack()}: }{Set back parts of frames.}
#' \item{\code{getWindowTitle()}: }{Get the title of the window.}
#' \item{\code{getHelp()}: }{Get the title of the help document.}
#' \item{\code{getParms()}: }{Get graphics settings parameters.}
#' \item{\code{checkTheme(index)}: }{Check themes.}
#' \item{\code{checkVariable(var)}: }{Check a variable length.}
#' \item{\code{checkError(parms)}: }{Check errors.}
#' \item{\code{setDataframe(parms)}: }{Set data frames.}
#' \item{\code{getGgplot(parms)}: }{Get \code{ggplot}.}
#' \item{\code{getGeom(parms)}: }{Get \code{geom}.}
#' \item{\code{getScale(parms)}: }{Get \code{scale}.}
#' \item{\code{getCoord(parms)}: }{Get \code{coord}.}
#' \item{\code{getFacet(parms)}: }{Get \code{facet}.}
#' \item{\code{getXlab(parms)}: }{Get \code{xlab}.}
#' \item{\code{getYlab(parms)}: }{Get \code{ylab}.}
#' \item{\code{getZlab(parms)}: }{Get \code{zlab}.}
#' \item{\code{getMain(parms)}: }{Get the main label.}
#' \item{\code{getTheme(parms)}: }{Get \code{theme}.}
#' \item{\code{getOpts(parms)}: }{Get other \code{opts}.}
#' \item{\code{getPlot(parms)}: }{Get the plot object.}
#' \item{\code{getMessage()}: }{Get the plot error message.}
#' \item{\code{commandDoIt(command)}: }{An wrapper function for command execution.}
#' }
#' @family plot
#'
#' @export plot_base
#' @name plot_base-class
#' @aliases plot_base
#' @rdname plot-plot_base
#' @docType class
#' @importFrom ggthemes theme_tufte theme_economist theme_solarized theme_stata
#' @importFrom ggthemes theme_excel theme_igray theme_few theme_calc theme_fivethirtyeight
#' @importFrom ggthemes theme_gdocs theme_hc theme_pander
#' @importFrom grid unit
#' @keywords hplot
plot_base <- setRefClass(

  Class = "plot_base",

  fields = c("top", "alternateFrame", "rmlist", "codes", "mode"),

  methods = list(

    #' Plot Windows
    plotWindow = function() {

      # note: The initializeDialog() generates "top"
      initializeDialog(window = topwindow, title = getWindowTitle())
      top            <<- topwindow
      alternateFrame <<- tkframe(top)

      setFront()

      # OK
      onOK <- function() {

        # doItAndPrint mode
        mode <<- 1
        parms <- getParms()

        closeDialog()

        errorCode <- checkError(parms)
        if (errorCode == TRUE) {
          removeRmlist()
          return()
        } else if (errorCode == FALSE) {

          logger("require(\"ggplot2\")")

          setDataframe(parms)
          
          .plot <- getPlot(parms)
          logger("print(.plot)")
          response <- tryCatch({
              print(.plot)
              ""
            }, error = function(ex) {
              tclvalue(RcmdrTkmessageBox(
                message = getMessage(),
                title   = gettextKmg2("Error"),
                icon    = "error",
                type    = "ok",
                default = "ok"
              ))
            }
          )
          if (response == "ok") {
            removeRmlist()
            return()
          }
          if (parms$save == "1") savePlot(.plot)
        }
        
        removeRmlist()

        activateMenus()
        tkfocus(CommanderWindow())

      }

      setBack()
      

      # note: The OKCancelHelp() generates "buttonsFrame"
      OKCancelHelp(window = top, helpSubject = getHelp())
      
      # Preview
      onPreview <- function() {
        
        # justDoIt mode
        mode <<- 0

        parms <- getParms()

        errorCode <- checkError(parms)
        if (errorCode == TRUE) {
          removeRmlist()
          return()
        } else if (errorCode == FALSE) {
          
          setDataframe(parms)
          
          .plot <- getPlot(parms)
          response <- tryCatch({
            print(.plot)
            ""
          }, error = function(ex) {
            tclvalue(RcmdrTkmessageBox(
              message = getMessage(),
              title   = gettextKmg2("Error"),
              icon    = "error",
              type    = "ok",
              default = "ok"
            ))
          }
          )
          if (response == "ok") {
            removeRmlist()
            return()
          }
        }
        removeRmlist()
        
      }
      previewButton <- buttonRcmdr(
        rightButtonsBox, text = gettextKmg2("Preview"), foreground = "yellow",
        width = nchar(gettextKmg2("Preview")), command = onPreview, 
        image = "::image::applyIcon", compound = "left"
        )

      tkgrid(previewButton, row = 0, column = 3, sticky = "nw")
      tkgrid.configure(previewButton, padx = c(5, 0))
      tkgrid(buttonsFrame, sticky = "nw")
      dialogSuffix()

      return()

    },

    #' Save Plot
    savePlot = function(plot) {

      plotName <- deparse(substitute(plot))
      if (.Platform$OS.type == "windows") {
        file <- tclvalue(tkgetSaveFile(
          filetypes = paste("{{All Files} {*}}",
                            "{{pdf (Portable Document Format)} {.pdf}}",
                            "{{jpg (Joint Photographic Experts Group)} {.jpg}}",
                            "{{tiff (Tagged Image File Format)} {.tiff}}",
                            "{{bmp (Bitmap Image)} {.bmp}}",
                            "{{svg (Scalable Vector Graphics)} {.svg}}",
                            "{{png (Portable Network Graphics)} {.png}}"),
          defaultextension = ".png", initialfile = "Rplots.png"
        ))
      } else {
        file <- tclvalue(tkgetSaveFile(
          filetypes = paste("{{png (Portable Network Graphics)} {.png}}",
                            "{{pdf (Portable Document Format)} {.pdf}}",
                            "{{jpg (Joint Photographic Experts Group)} {.jpg}}",
                            "{{tiff (Tagged Image File Format)} {.tiff}}",
                            "{{bmp (Bitmap Image)} {.bmp}}",
                            "{{svg (Scalable Vector Graphics)} {.svg}}"),
          initialfile = "Rplots"
        ))
      }
      if (file == "") return()

      if (class(.self)[1] == "gkm") {
        command <- paste0("RcmdrPlugin.KMggplot2::ggsaveKmg2(filename = \"", file, "\", plot = ", plotName, ")")
      } else {
        command <- paste0("ggsave(filename = \"", file, "\", plot = ", plotName, ")")
      }
      doItAndPrint(command)

      return()

    },

    #' Register \code{rm()} List
    registRmlist = function(object) {
      
      txtObjects <- deparse(substitute(object))
      if (class(rmlist) == "uninitializedField") {
        rmlist <<- list(txtObjects)
      } else {
        rmlist <<- c(rmlist, txtObjects)
      }
      return()

    },

    #' Remove \code{rm()} List
    removeRmlist = function() {

      if (class(rmlist) != "uninitializedField") {
        command <- do.call(paste, c(rmlist, sep=", "))
        command <- paste0("rm(", command, ")")
        commandDoIt(command)
        rmlist <<- list()
      }
      return()

    },

    #' Set Front
    setFront = function() {

      return()

    },

    #' Set Back
    setBack = function() {

      return()

    },

    #' Get Window Title
    getWindowTitle = function() {

      "plot_base"

    },

    #' Get Help
    getHelp = function() {
      
      "plot_base"
      
    },

    #' Get Parameters
    getParms = function() {

      x      <- "x"
      y      <- "y"
      z      <- "z"
      s      <- "s"
      t      <- "t"

      x      <- checkVariable(x)
      y      <- checkVariable(y)
      z      <- checkVariable(z)
      s      <- checkVariable(s)
      t      <- checkVariable(t)

      xlab   <- ""
      xauto  <- x
      ylab   <- ""
      yauto  <- y
      zlab   <- ""
      main   <- ""

      size   <- "16"
      family <- "0"
      colour <- "0"
      save   <- "0"
      theme  <- "0"

      options(
        kmg2FontSize   = size,
        kmg2FontFamily = family,
        kmg2ColourSet  = colour,
        kmg2SaveGraph  = save,
        kmg2Theme      = theme
      )

      list(
        x = x, y = y, z = z, s = s, t = t,
        xlab = xlab, xauto = xauto, ylab = ylab, yauto = yauto, zlab = zlab, main = main,
        size = size, family = family, colour = colour, save = save, theme = theme
      )

    },

    #' Check themes.
    checkTheme = function(index) {

      if (index == "theme_bw") {
        theme <- "theme_bw"
      } else if (index == "theme_simple") {
        theme <- "RcmdrPlugin.KMggplot2::theme_simple"
      } else if (index == "theme_classic") {
        theme <- "theme_classic"
      } else if (index == "theme_gray") {
        theme <- "theme_gray"
      } else if (index == "theme_minimal") {
        theme <- "theme_minimal"
      } else if (index == "theme_linedraw") {
        theme <- "theme_linedraw"
      } else if (index == "theme_light") {
        theme <- "theme_light"
      } else if (index == "theme_dark") {
        theme <- "theme_dark"
      } else if (index == "theme_base") {
        theme <- "ggthemes::theme_base"
        commandDoIt("ggthemes_data <- ggthemes::ggthemes_data")
        registRmlist(ggthemes_data)
      } else if (index == "theme_calc") {
        theme <- "ggthemes::theme_calc"
        commandDoIt("ggthemes_data <- ggthemes::ggthemes_data")
        registRmlist(ggthemes_data)
      } else if (index == "theme_economist") {
        theme <- "ggthemes::theme_economist"
        commandDoIt("ggthemes_data <- ggthemes::ggthemes_data")
        registRmlist(ggthemes_data)
      } else if (index == "theme_excel") {
        theme <- "ggthemes::theme_excel"
      } else if (index == "theme_few") {
        theme <- "ggthemes::theme_few"
        commandDoIt("ggthemes_data <- ggthemes::ggthemes_data")
        registRmlist(ggthemes_data)
      } else if (index == "theme_fivethirtyeight") {
        theme <- "ggthemes::theme_fivethirtyeight"
        commandDoIt("ggthemes_data <- ggthemes::ggthemes_data")
        registRmlist(ggthemes_data)
      } else if (index == "theme_gdocs") {
        theme <- "ggthemes::theme_gdocs"
      } else if (index == "theme_hc") {
        theme <- "ggthemes::theme_hc"
        commandDoIt("ggthemes_data <- ggthemes::ggthemes_data")
        registRmlist(ggthemes_data)
      } else if (index == "theme_par") {
        theme <- "ggthemes::theme_par"
        commandDoIt("ggthemes_data <- ggthemes::ggthemes_data")
        registRmlist(ggthemes_data)
      } else if (index == "theme_pander") {
        theme <- "ggthemes::theme_pander"
      } else if (index == "theme_solarized") {
        theme <- "ggthemes::theme_solarized"
        commandDoIt("ggthemes_data <- ggthemes::ggthemes_data")
        registRmlist(ggthemes_data)
      } else if (index == "theme_stata") {
        theme <- "ggthemes::theme_stata"
        commandDoIt("ggthemes_data <- ggthemes::ggthemes_data")
        registRmlist(ggthemes_data)
      } else if (index == "theme_tufte") {
        theme <- "ggthemes::theme_tufte"
      } else if (index == "theme_wsj2") {
        theme <- "RcmdrPlugin.KMggplot2::theme_wsj2"
        commandDoIt("ggthemes_data <- ggthemes::ggthemes_data")
        registRmlist(ggthemes_data)
      } else if (index == "theme_igray") {
        theme <- "ggthemes::theme_igray"
      } else {
        theme <- "theme_bw"
      }
      theme
      
    },

    #' Check variable length.
    checkVariable = function(var) {
      
      if (length(var) > 1) {
        var <- var[1]
      }
      var
      
    },

    #' Check Error
    checkError = function(parms) {

      errorCode <- FALSE

    },

    #' Set \code{data.frame}
    setDataframe = function(parms) {
      
      var <- list()
      if (length(parms$x) != 0) {
        var <- c(var, paste0("x = ", ActiveDataSet(), "$", parms$x))
      }
      if (length(parms$y) != 0) {
        var <- c(var, paste0("y = ", ActiveDataSet(), "$", parms$y))
      }
      if (length(parms$z) != 0) {
        var <- c(var, paste0("z = ", ActiveDataSet(), "$", parms$z))
      }
      if (length(parms$s) != 0) {
        var <- c(var, paste0("s = ", ActiveDataSet(), "$", parms$s))
      }
      if (length(parms$t) != 0) {
        var <- c(var, paste0("t = ", ActiveDataSet(), "$", parms$t))
      }
      command <- do.call(paste, c(var, list(sep = ", ")))
      command <- paste0(".df <- data.frame(", command, ")")

      commandDoIt(command)
      registRmlist(.df)

    },
    
    #' Get Ggplot
    getGgplot = function(parms) {
      
      "ggplot(data.frame(1), aes(x = 1, y = 1)) + "
      
    },
    
    #' Get Geom
    getGeom = function(parms) {
      
      "geom_point() + "
      
    },
    
    #' Get Scale
    getScale = function(parms) {
      
      "scale_y_continuous(expand = c(0.01, 0)) + "
      
    },
    
    #' Get Coord
    getCoord = function(parms) {
      
      ""
      
    },

    #' Get Facet
    getFacet = function(parms) {

      if (length(parms$s) != 0 && length(parms$t) != 0) {
        facet <- "facet_grid(s ~ t) + "
      } else if (length(parms$s) != 0) {
        facet <- "facet_wrap( ~ s) + "
      } else if (length(parms$t) != 0) {
        facet <- "facet_wrap( ~ t) + "
      } else {
        facet <- ""
      }
      facet

    },
    
    #' Get Xlab
    getXlab = function(parms) {

      if (nchar(parms$xlab) == 0) {
        xlab <- "xlab(NULL) + "
      } else if (parms$xlab == "<auto>") {
        xlab <- paste0("xlab(\"", parms$xauto, "\") + ")
      } else {
        xlab <- paste0("xlab(\"", parms$xlab, "\") + ")
      }
      xlab

    },
    
    #' Get Ylab
    getYlab = function(parms) {

      if (nchar(parms$ylab) == 0) {
        ylab <- "ylab(NULL) + "
      } else if (parms$ylab == "<auto>") {
        ylab <- paste0("ylab(\"", parms$yauto, "\") + ")
      } else {
        ylab <- paste0("ylab(\"", parms$ylab, "\") + ")
      }
      ylab

    },
    
    #' Get Zlab
    getZlab = function(parms) {

      ""

    },

    #' Get Main
    getMain = function(parms) {

      if (nchar(parms$main) == 0) {
        main <- ""
      } else {
        main <- paste0("labs(title = \"", parms$main, "\") + ")
      }
      main

    },
    
    #' Get Theme
    getTheme = function(parms) {

      paste0(parms$theme, "(base_size = ", parms$size, ", base_family = \"", parms$family, "\")")

    },
    
    #' Get Opts
    getOpts = function(parms) {
      
      opts <- list()
      if (length(parms$s) != 0 || length(parms$t) != 0) {
        opts <- c(opts, "panel.margin = grid::unit(0.3, \"lines\")")
      }

      if (length(opts) != 0) {
        opts <- do.call(paste, c(opts, list(sep = ", ")))
        opts <- paste0(" + theme(", opts, ")")
      } else {
        opts <- ""
      }
      opts
      
    },
    
    #' Get Plot
    getPlot = function(parms) {

      gg    <- getGgplot(parms)
      geom  <- getGeom(parms)
      scale <- getScale(parms)
      coord <- getCoord(parms)
      facet <- getFacet(parms)
      xlab  <- getXlab(parms)
      ylab  <- getYlab(parms)
      zlab  <- getZlab(parms)
      main  <- getMain(parms)
      theme <- getTheme(parms)
      opts  <- getOpts(parms)

      command <- paste0(
        ".plot <- ",
        gg, geom, scale, coord, facet,
        xlab, ylab, zlab, main, theme, opts
      )
      commandDoIt(command)
      registRmlist(.plot)
      return(.plot)

    },

    #' Get Plot Error Message
    getMessage = function() {

      gettextKmg2("Plot failed.  Please check the data and variables, or try other options.")

    },
    
    #' An wrapper function for command execution
    commandDoIt = function(command, log = TRUE) {
      
      if (mode == 1) {
        doItAndPrint(command, log = log)
      } else {
        justDoIt(command)
      }
      NULL
      
    }

  )
)



#' Create Plot Window
#'
#' \code{plotWindow} method creates the window that make plots.
#'
#' @usage \S4method{plotWindow}{plot_base}()
#' @family plot
#'
#' @name plotWindow,plot_base-method
#' @rdname plot-plot_base-plotWindow
#' @docType methods
#' @keywords hplot
NULL



#' Save Plot
#'
#' \code{savePlot} method saves the plot.
#'
#' @usage \S4method{savePlot}{plot_base}(plot)
#' @param plot \code{ggplot} or \code{recordedplot} class object; the plot to save.
#' @family plot
#'
#' @name savePlot,plot_base-method
#' @rdname plot-plot_base-savePlot
#' @docType methods
#' @keywords hplot
NULL



#' Register Deletable Temporary Objects
#'
#' \code{registRmlist} method registers deletable temporary objects.
#'
#' @usage \S4method{registRmlist}{plot_base}(object)
#' @param object Object; the deletable temporary object to register.
#' @family plot
#'
#' @name registRmlist,plot_base-method
#' @rdname plot-plot_base-registRmlist
#' @docType methods
#' @keywords hplot
NULL



#' Remove Registered Temporary Objects
#'
#' \code{removeRmlist} method removes registered temporary objects.
#'
#' @usage \S4method{removeRmlist}{plot_base}()
#' @family plot
#'
#' @name removeRmlist,plot_base-method
#' @rdname plot-plot_base-removeRmlist
#' @docType methods
#' @keywords hplot
NULL



#' Set Front Parts
#'
#' \code{setFront} method sets front parts of frames.
#'
#' @usage \S4method{setFront}{plot_base}()
#' @family plot
#'
#' @name setFront,plot_base-method
#' @rdname plot-plot_base-setFront
#' @docType methods
#' @keywords hplot
NULL



#' Set Back Parts
#'
#' \code{setBack} method sets back parts of frames.
#'
#' @usage \S4method{setBack}{plot_base}()
#' @family plot
#'
#' @name setBack,plot_base-method
#' @rdname plot-plot_base-setBack
#' @docType methods
#' @keywords hplot
NULL



#' Get Title
#'
#' \code{getWindowTitle} method gets the title of the window.
#'
#' @usage \S4method{getWindowTitle}{plot_base}()
#' @return
#' \describe{
#' \item{title}{Character; the title of the window.}
#' }
#' @family plot
#'
#' @name getWindowTitle,plot_base-method
#' @rdname plot-plot_base-getWindowTitle
#' @docType methods
#' @keywords hplot
NULL



#' Get Help
#'
#' \code{getHelp} method gets the title of the help document.
#'
#' @usage \S4method{getHelp}{plot_base}()
#' @return
#' \describe{
#' \item{help}{Character; the title of the help document.}
#' }
#' @family plot
#'
#' @name getHelp,plot_base-method
#' @rdname plot-plot_base-getHelp
#' @docType methods
#' @keywords hplot
NULL



#' Get Graphics Parameters
#'
#' \code{getParms} method gets graphics settings parameters.
#'
#' @usage \S4method{getParms}{plot_base}()
#' @return
#' \describe{
#' \item{parms}{List of objects; graphics settings parameters.}
#' }
#' @family plot
#'
#' @name getParms,plot_base-method
#' @rdname plot-plot_base-getParms
#' @docType methods
#' @keywords hplot
NULL



#' Check Themes
#'
#' \code{checkTheme} method checks a theme index.
#'
#' @usage \S4method{checkTheme}{plot_base}(index)
#' @param index Character; a theme index.
#' @return
#' \describe{
#' \item{theme}{Character; the theme name.}
#' }
#' @family plot
#'
#' @name checkTheme,plot_base-method
#' @rdname plot-plot_base-checkTheme
#' @docType methods
#' @keywords hplot
NULL



#' Check Variable Length
#'
#' \code{checkVariable} method checks a variable length.
#'
#' @usage \S4method{checkVariable}{plot_base}(var)
#' @param var Vector of Character; variable names.
#' @return
#' \describe{
#' \item{var}{Character; the fixed variable name.}
#' }
#' @family plot
#'
#' @name checkVariable,plot_base-method
#' @rdname plot-plot_base-checkVariable
#' @docType methods
#' @keywords hplot
NULL



#' Check Errors
#'
#' \code{checkError} method checks errors.
#'
#' @usage \S4method{checkError}{plot_base}(parms)
#' @param parms List of objects; graphics settings parameters.
#' @return
#' \describe{
#' \item{errorCode}{Boolean; the error code.}
#' }
#' @family plot
#'
#' @name checkError,plot_base-method
#' @rdname plot-plot_base-checkError
#' @docType methods
#' @keywords hplot
NULL



#' Set Data Frames
#'
#' \code{setDataframe} method sets data frames.
#'
#' @usage \S4method{setDataframe}{plot_base}(parms)
#' @param parms List of objects; graphics settings parameters.
#' @family plot
#'
#' @name setDataframe,plot_base-method
#' @rdname plot-plot_base-setDataframe
#' @docType methods
#' @keywords hplot
NULL



#' Get \code{ggplot}
#'
#' \code{getGgplot} method gets \code{ggplot}.
#'
#' @usage \S4method{getGgplot}{plot_base}(parms)
#' @param parms List of objects; graphics settings parameters.
#' @return
#' \describe{
#' \item{ggplot}{Character; a code snippet for ggplot functions.}
#' }
#' @family plot
#'
#' @name getGgplot,plot_base-method
#' @rdname plot-plot_base-getGgplot
#' @docType methods
#' @keywords hplot
NULL



#' Get \code{geom}
#'
#' \code{getGeom} method gets \code{geom}.
#'
#' @usage \S4method{getGeom}{plot_base}(parms)
#' @param parms List of objects; graphics settings parameters.
#' @return
#' \describe{
#' \item{geom}{Character; a code snippet for geometric objects.}
#' }
#' @family plot
#'
#' @name getGeom,plot_base-method
#' @rdname plot-plot_base-getGeom
#' @docType methods
#' @keywords hplot
NULL



#' Get \code{scale}
#'
#' \code{getScale} method gets \code{scale}.
#'
#' @usage \S4method{getScale}{plot_base}(parms)
#' @param parms List of objects; graphics settings parameters.
#' @return
#' \describe{
#' \item{scale}{Character; a code snippet for scales.}
#' }
#' @family plot
#'
#' @name getScale,plot_base-method
#' @rdname plot-plot_base-getScale
#' @docType methods
#' @keywords hplot
NULL



#' Get \code{coord}
#'
#' \code{getCoord} method gets \code{coord}.
#'
#' @usage \S4method{getCoord}{plot_base}(parms)
#' @param parms List of objects; graphics settings parameters.
#' @return
#' \describe{
#' \item{coord}{Character; a code snippet for coords.}
#' }
#' @family plot
#'
#' @name getCoord,plot_base-method
#' @rdname plot-plot_base-getCoord
#' @docType methods
#' @keywords hplot
NULL



#' Get \code{facet}
#'
#' \code{getFacet} method gets \code{facet}.
#'
#' @usage \S4method{getFacet}{plot_base}(parms)
#' @param parms List of objects; graphics settings parameters.
#' @return
#' \describe{
#' \item{facet}{Character; a code snippet for the facet.}
#' }
#' @family plot
#'
#' @name getFacet,plot_base-method
#' @rdname plot-plot_base-getFacet
#' @docType methods
#' @keywords hplot
NULL



#' Get \code{xlab}
#'
#' \code{getXlab} method gets the x-axis label.
#'
#' @usage \S4method{getXlab}{plot_base}(parms)
#' @param parms List of objects; graphics settings parameters.
#' @return
#' \describe{
#' \item{xlab}{Character; a code snippet for the x-axis label.}
#' }
#' @family plot
#'
#' @name getXlab,plot_base-method
#' @rdname plot-plot_base-getXlab
#' @docType methods
#' @keywords hplot
NULL



#' Get \code{ylab}
#'
#' \code{getYlab} method gets the y-axis label.
#'
#' @usage \S4method{getYlab}{plot_base}(parms)
#' @param parms List of objects; graphics settings parameters.
#' @return
#' \describe{
#' \item{ylab}{Character; a code snippet for the y-axis label.}
#' }
#' @family plot
#'
#' @name getYlab,plot_base-method
#' @rdname plot-plot_base-getYlab
#' @docType methods
#' @keywords hplot
NULL



#' Get \code{zlab}
#'
#' \code{getZlab} method gets the label of the stratum.
#'
#' @usage \S4method{getZlab}{plot_base}(parms)
#' @param parms List of objects; graphics settings parameters.
#' @return
#' \describe{
#' \item{zlab}{Character; a code snippet for the label of the stratum.}
#' }
#' @family plot
#'
#' @name getZlab,plot_base-method
#' @rdname plot-plot_base-getZlab
#' @docType methods
#' @keywords hplot
NULL



#' Get Main Label
#'
#' \code{getMain} method gets the main label (graph title).
#'
#' @usage \S4method{getMain}{plot_base}(parms)
#' @param parms List of objects; graphics settings parameters.
#' @return
#' \describe{
#' \item{main}{Character; a code snippet for the main label.}
#' }
#' @family plot
#'
#' @name getMain,plot_base-method
#' @rdname plot-plot_base-getMain
#' @docType methods
#' @keywords hplot
NULL



#' Get \code{theme}
#'
#' \code{getTheme} method gets the graph \code{theme}.
#'
#' @usage \S4method{getTheme}{plot_base}(parms)
#' @param parms List of objects; graphics settings parameters.
#' @return
#' \describe{
#' \item{theme}{Character; a code snippet for the \code{theme}.}
#' }
#' @family plot
#'
#' @name getTheme,plot_base-method
#' @rdname plot-plot_base-getTheme
#' @docType methods
#' @keywords hplot
NULL



#' Get Other \code{opts}
#'
#' \code{getOpts} method gets other graph options.
#'
#' @usage \S4method{getOpts}{plot_base}(parms)
#' @param parms List of objects; graphics settings parameters.
#' @return
#' \describe{
#' \item{opts}{Character; a code snippet for other graph options.}
#' }
#' @family plot
#'
#' @name getOpts,plot_base-method
#' @rdname plot-plot_base-getOpts
#' @docType methods
#' @keywords hplot
NULL



#' Get Plot
#'
#' \code{getPlot} method gets the plot object.
#'
#' @usage \S4method{getPlot}{plot_base}(parms)
#' @param parms List of objects; graphics settings parameters.
#' @return
#' \describe{
#' \item{plot}{\code{ggplot2} class object; the plot.}
#' }
#' @family plot
#'
#' @name getPlot,plot_base-method
#' @rdname plot-plot_base-getPlot
#' @docType methods
#' @keywords hplot
NULL



#' Get Plot Error Message
#'
#' \code{getMessage} method gets the plot error message.
#'
#' @usage \S4method{getMessage}{plot_base}()
#' @return
#' \describe{
#' \item{message}{Character; the plot error message.}
#' }
#' @family plot
#'
#' @name getMessage,plot_base-method
#' @rdname plot-plot_base-getMessage
#' @docType methods
#' @keywords hplot
NULL



#' An Wrapper Function for Command Execution
#'
#' \code{commandDoIt} method gets the plot error message.
#'
#' @usage \S4method{commandDoIt}{plot_base}(command)
#' @param command String; command codes.
#' @family plot
#'
#' @name commandDoIt,plot_base-method
#' @rdname plot-plot_base-commandDoIt
#' @docType methods
#' @keywords hplot
NULL
