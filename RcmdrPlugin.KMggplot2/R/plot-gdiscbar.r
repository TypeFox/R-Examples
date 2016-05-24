#' Discrete Bar Charts Subclass
#'
#' \code{gdiscbar} class is a subclass for discrete bar charts.
#'
#' This class is a subclass which show dialog boxes of discrete bar charts for graphics editing.
#'
#' @section Fields:
#' \describe{
#' \item{\code{top}: }{\code{tkwin} class object; parent of widget window.} 
#' \item{\code{alternateFrame}: }{\code{tkwin} class object; a special frame for some GUI parts.} 
#' \item{\code{vbbox1}: }{\code{variableboxes} class object; the frame to select variables.} 
#' \item{\code{vbbox2}: }{\code{variableboxes} class object; the frame to select facet variables.} 
#' \item{\code{lbbox1}: }{\code{textfields} class object; the frame to set axis labels and the main title.} 
#' \item{\code{rbbox1}: }{\code{radioboxes} class object; the frame to set the axis scaling.} 
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
#' @name gdiscbar-class
#' @aliases gdiscbar
#' @rdname plot-gdiscbar
#' @docType class
#' @keywords hplot
#' @importFrom scales percent_format
#' @export gdiscbar
gdiscbar <- setRefClass(

  Class = "gdiscbar",

  fields = c("vbbox1", "vbbox2", "lbbox1", "rbbox1", "tbbox1"),

  contains = c("plot_base"),

  methods = list(

    setFront = function() {

      vbbox1 <<- variableboxes$new()
      vbbox1$front(
        top       = top, 
        types     = list(Factors(), Factors()),
        titles    = list(
          gettextKmg2("X variable (pick one)"),
          gettextKmg2("Stratum variable")
        ),
        initialSelection = list(0, FALSE)
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
        initValues = list("<auto>", "<auto>", "<auto>", ""),
        titles     = list(
          gettextKmg2("Horizontal axis label"),
          gettextKmg2("Vertical axis label"),
          gettextKmg2("Legend label"),
          gettextKmg2("Title")
        )
      )

      rbbox1 <<- radioboxes$new()
      rbbox1$front(
        top    = top,
        labels = list(
          gettextKmg2("Percentages"),
          gettextKmg2("Frequency counts")
        ),
        title  = gettextKmg2("Axis scaling")
      )

      tbbox1 <<- toolbox$new()
      tbbox1$front(top)

    },

    setBack = function() {

      vbbox1$back()
      vbbox2$back()
      lbbox1$back()
      rbbox1$back()
      tbbox1$back()

    },

    getWindowTitle = function() {
      
      gettextKmg2("Bar chart for discrete variables")
      
    },
    
    getHelp = function() {
      
      "geom_bar"
      
    },

    getParms = function() {

      x      <- getSelection(vbbox1$variable[[1]])
      y      <- character(0)
      z      <- getSelection(vbbox1$variable[[2]])

      s      <- getSelection(vbbox2$variable[[1]])
      t      <- getSelection(vbbox2$variable[[2]])

      x      <- checkVariable(x)
      z      <- checkVariable(z)
      s      <- checkVariable(s)
      t      <- checkVariable(t)

      xlab   <- tclvalue(lbbox1$fields[[1]]$value)
      xauto  <- x
      ylab   <- tclvalue(lbbox1$fields[[2]]$value)
      # yauto  <- y
      zlab   <- tclvalue(lbbox1$fields[[3]]$value)
      main   <- tclvalue(lbbox1$fields[[4]]$value)

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
      
      axisScaling  <- tclvalue(rbbox1$value)
      if (axisScaling == "1") {
        yauto  <- gettextKmg2("Percent")
      } else {
        yauto  <- gettextKmg2("Count")
      }

      list(
        x = x, y = y, z = z, s = s, t = t,
        xlab = xlab, xauto = xauto, ylab = ylab, yauto = yauto, zlab = zlab, main = main,
        size = size, family = family, colour = colour, save = save, theme = theme,
        axisScaling = axisScaling
      )

    },

    checkError = function(parms) {

      if (length(parms$x) == 0) {
        errorCondition(
          recall  = windowDiscretebar,
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
      table <- list()
      if (length(parms$x) != 0) {
        var <- c(var, paste0("x = ", ActiveDataSet(), "$", parms$x))
        table <- c(table, "x")
      }
      if (length(parms$z) != 0) {
        var <- c(var, paste0("z = ", ActiveDataSet(), "$", parms$z))
        table <- c(table, "z")
      }
      if (length(parms$s) != 0) {
        var <- c(var, paste0("s = ", ActiveDataSet(), "$", parms$s))
        table <- c(table, "s")
      }
      if (length(parms$t) != 0) {
        var <- c(var, paste0("t = ", ActiveDataSet(), "$", parms$t))
        table <- c(table, "t")
      }
      command <- do.call(paste, c(var, list(sep = ", ")))
      command <- paste0(".df <- data.frame(", command, ")")


      if (length(parms$s) != 0 && length(parms$t) != 0) {
        margin <- paste0(length(table) - 1, ":", length(table))
      } else if (length(parms$s) != 0 || length(parms$t) != 0) {
        margin <- paste0(length(table))
      } else {
        margin <- "NULL"
      }
      table <- do.call(paste, c(table, list(sep = ", ")))

      if (parms$axisScaling == "1") {
        table <- paste0("prop.table(table(", table, "), margin = ", margin, ")")
      } else {
        table <- paste0("table(", table, ")")
      }

      command <- paste0(
        command,
        "\n",
        ".df <- as.data.frame(with(.df, ", table, "))"
      )

      commandDoIt(command)
      registRmlist(.df)

    },

    getGgplot = function(parms) {

      if (length(parms$z) != 0) {
        ggplot <- "ggplot(data = .df, aes(x = x, y = Freq, fill = z)) + "
      } else {
        ggplot <- "ggplot(data = .df, aes(x = x, y = Freq)) + "
      }
      ggplot

    },

    getGeom = function(parms) {

      if (length(parms$z) != 0) {
        if (parms$axisScaling == "1") {
          geom <- "geom_bar(width = 0.9, position = \"fill\", stat = \"identity\") + "
        } else {
          geom <- "geom_bar(width = 0.9, position = \"stack\", stat = \"identity\") + "
        }
      } else {
        geom <- "geom_bar(width = 0.9, stat = \"identity\") + "
      }
      geom

    },

    getScale = function(parms) {
      
      if (length(parms$z) != 0) {
        scale <- paste0("scale_fill_brewer(palette = \"", parms$colour, "\") + ")
      } else {
        scale <- ""
      }

      if (parms$axisScaling == "1") {
        scale <- paste0(
          scale,
          "scale_y_continuous(expand = c(0.01, 0), labels = scales::percent_format()) + "
        )
      } else {
        scale <- paste0(
          scale,
          "scale_y_continuous(expand = c(0.01, 0)) + "
        )
      }
      scale

    },

    getZlab = function(parms) {

      if (length(parms$z) == 0) {
        zlab <- ""
      } else if (nchar(parms$zlab) == 0) {
        zlab <- ""
      } else if (parms$zlab == "<auto>") {
        zlab <- paste0("labs(fill = \"", parms$z, "\") + ")
      } else {
        zlab <- paste0("labs(fill = \"", parms$zlab, "\") + ")
      }
      zlab

    },

    getOpts = function(parms) {

      opts <- list()
      if (length(parms$s) != 0 || length(parms$t) != 0) {
        opts <- c(opts, "panel.margin = unit(0.3, \"lines\")")
      }

      if (length(parms$z) != 0 && nchar(parms$zlab) == 0) {
        opts <- c(opts, "legend.position = \"right\"", "legend.title = element_blank()")
      } else if (length(parms$z) != 0 && nchar(parms$zlab) != 0) {
        opts <- c(opts, "legend.position = \"right\"")
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



#' Wrapper Function of Discrete Bar Charts Subclass
#'
#' \code{windowDiscretebar} function is a wrapper function of \code{gdiscbar} class for the R-commander menu bar.
#'
#' @rdname plot-gdiscbar-windowDiscretebar
#' @keywords hplot
#' @export
windowDiscretebar <- function() {

  Discretebar <- RcmdrPlugin.KMggplot2::gdiscbar$new()
  Discretebar$plotWindow()

}
