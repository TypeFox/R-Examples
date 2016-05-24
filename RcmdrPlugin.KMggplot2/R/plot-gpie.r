#' Pie Charts Subclass
#'
#' \code{gpie} class is a subclass for pie charts.
#'
#' This class is a subclass which show dialog boxes of pie charts for graphics editing.
#'
#' @section Fields:
#' \describe{
#' \item{\code{top}: }{\code{tkwin} class object; parent of widget window.} 
#' \item{\code{alternateFrame}: }{\code{tkwin} class object; a special frame for some GUI parts.} 
#' \item{\code{vbbox1}: }{\code{variableboxes} class object; the frame to select variables.} 
#' \item{\code{vbbox2}: }{\code{variableboxes} class object; the frame to select facet variables.} 
#' \item{\code{lbbox1}: }{\code{textfields} class object; the frame to set axis labels and the main title.} 
#' \item{\code{tbbox1}: }{\code{toolbox} class object; the frame to set the font, the colour set, other option, and the theme.}
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
#' @name gpie-class
#' @aliases gpie
#' @rdname plot-gpie
#' @docType class
#' @keywords hplot
#' @importFrom scales percent_format
#' @export gpie
gpie <- setRefClass(

  Class = "gpie",

  fields = c("vbbox1", "vbbox2", "lbbox1", "tbbox1"),

  contains = c("plot_base"),

  methods = list(

    setFront = function() {

      vbbox1 <<- variableboxes$new()
      vbbox1$front(
        top       = top, 
        types     = list(Factors()),
        titles    = list(
          gettextKmg2("Variable (pick one)")
        ),
        initialSelection = list(0)
      )

      vbbox2 <<- variableboxes$new()
      vbbox2$front(
        top       = top, 
        types     = list(Factors(), Factors()),
        titles    = list(
          gettextKmg2("Facet variable in rows"),
          gettextKmg2("Facet variable in cols")
        )
      )

      lbbox1 <<- textfields$new()
      lbbox1$front(
        top        = top,
        initValues = list("<auto>", ""),
        titles     = list(
          gettextKmg2("Legend label"),
          gettextKmg2("Title")
        )
      )

      tbbox1 <<- toolbox$new()
      tbbox1$front(top)

    },

    setBack = function() {

      vbbox1$back()
      vbbox2$back()
      lbbox1$back()
      tbbox1$back()

    },

    getWindowTitle = function() {
      
      gettextKmg2("Pie chart")
      
    },
    
    getHelp = function() {
      
      "coord_polar"
      
    },

    getParms = function() {

      x      <- character(0)
      y      <- getSelection(vbbox1$variable[[1]])
      z      <- character(0)

      s      <- getSelection(vbbox2$variable[[1]])
      t      <- getSelection(vbbox2$variable[[2]])

      y      <- checkVariable(y)
      s      <- checkVariable(s)
      t      <- checkVariable(t)

      xlab   <- ""
      xauto  <- x
      ylab   <- tclvalue(lbbox1$fields[[1]]$value)
      yauto  <- y
      zlab   <- ""
      main   <- tclvalue(lbbox1$fields[[2]]$value)

      size   <- tclvalue(tbbox1$size$value)
      family <- getSelection(tbbox1$family)
      colour <- getSelection(tbbox1$colour)
      save   <- tclvalue(tbbox1$goption$value[[1]])
      theme  <- checkTheme(getSelection(tbbox1$theme))
      
      options(
        kmg2FontSize   = tclvalue(tbbox1$size$value),
        kmg2FontFamily = seq_along(tbbox1$family$varlist)[tbbox1$family$varlist == getSelection(tbbox1$family)] - 1,
        kmg2ColourSet  = seq_along(tbbox1$colour$varlist)[tbbox1$colour$varlist == getSelection(tbbox1$colour)] - 1,
        kmg2SaveGraph  = tclvalue(tbbox1$goption$value[[1]]),
        kmg2Theme      = seq_along(tbbox1$theme$varlist)[tbbox1$theme$varlist == getSelection(tbbox1$theme)] - 1
      )
      
      list(
        x = x, y = y, z = z, s = s, t = t,
        xlab = xlab, ylab = ylab, zlab = zlab, main = main,
        size = size, family = family, colour = colour, save = save, theme = theme
      )

    },

    checkError = function(parms) {

      if (length(parms$y) == 0) {
        errorCondition(
          recall  = windowPie,
          message = gettextKmg2("Variable is not selected")
        )
        errorCode <- TRUE
      } else {
        errorCode <- FALSE
      }
      errorCode

    },

    setDataframe = function(parms) {

      var <- list()
      if (length(parms$y) != 0) {
        var <- c(var, paste0("y = ", ActiveDataSet(), "$", parms$y))
      }
      if (length(parms$s) != 0) {
        var <- c(var, paste0("s = ", ActiveDataSet(), "$", parms$s))
      }
      if (length(parms$t) != 0) {
        var <- c(var, paste0("t = ", ActiveDataSet(), "$", parms$t))
      }
      command <- do.call(paste, c(var, list(sep = ", ")))
      command <- paste0(".df <- data.frame(", command, ")")

      if (length(parms$s) != 0 && length(parms$t) != 0) {
        command <- paste0(
          command,
          "\n",
          ".df <- as.data.frame(with(.df, prop.table(table(y, s, t), margin = 2:3)))"
        )
      } else if (length(parms$s) != 0) {
        command <- paste0(
          command,
          "\n",
          ".df <- as.data.frame(with(.df, prop.table(table(y, s), margin = 2)))"
        )
      } else if (length(parms$t) != 0) {
        command <- paste0(
          command,
          "\n",
          ".df <- as.data.frame(with(.df, prop.table(table(y, t), margin = 2)))"
        )
      } else {
        command <- paste0(
          command,
          "\n",
          ".df <- as.data.frame(with(.df, prop.table(table(y))))"
        )
      }

      commandDoIt(command)
      registRmlist(.df)

    },

    getGgplot = function(parms) {

      "ggplot(data = .df, aes(x = factor(1), y = Freq, fill = y)) + "

    },

    getGeom = function(parms) {

      "geom_bar(width = 1, stat = \"identity\") + "

    },

    getScale = function(parms) {

      paste0(
        "scale_fill_brewer(palette = \"", parms$colour, "\") + ",
        "scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) + "
      )

    },

    getCoord = function(parms) {

      "coord_polar(theta = \"y\") + "

    },

    getYlab = function(parms) {

      if (nchar(parms$ylab) == 0) {
        ylab <- ""
      } else if (parms$ylab == "<auto>") {
        ylab <- paste0("labs(fill = \"", parms$y, "\") + ")
      } else {
        ylab <- paste0("labs(fill = \"", parms$ylab, "\") + ")
      }
      ylab

    },

    getOpts = function(parms) {

      opts <- list()
      if (length(parms$s) != 0 || length(parms$t) != 0) {
        opts <- c(opts, "panel.margin = unit(0.3, \"lines\")")
      }

      opts <- c(
        opts,
        "legend.position = \"right\"",
        "axis.title.x = element_blank()",
        "axis.title.y = element_blank()",
        "axis.text.y = element_blank()",
        "axis.ticks = element_blank()"
      )

      if (nchar(parms$ylab) == 0) {
        opts <- c(opts, "legend.title = element_blank()")
      }

      if (parms$theme == "RcmdrPlugin.KMggplot2::theme_simple") {
        opts <- c(opts, "panel.border = element_blank()")
      }

      if (length(opts) != 0) {
        opts <- do.call(paste, c(opts, list(sep = ", ")))
        opts <- paste0(" + theme(", opts, ")")
      } else {
        opts <- ""
      }
      opts

    }

  )
)



#' Wrapper Function of Pie Charts Subclass
#'
#' \code{windowPie} function is a wrapper function of \code{gpie} class for the R-commander menu bar.
#'
#' @rdname plot-gpie-windowPie
#' @keywords hplot
#' @export
windowPie <- function() {

  Pie <- RcmdrPlugin.KMggplot2::gpie$new()
  Pie$plotWindow()

}
