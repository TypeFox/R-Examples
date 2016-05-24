#' Line Chart Subclass
#'
#' \code{gline} class is a subclass for line charts.
#'
#' This class is a subclass which show dialog boxes of line charts for graphics editing.
#'
#' @section Fields:
#' \describe{
#' \item{\code{top}: }{\code{tkwin} class object; parent of widget window.}
#' \item{\code{alternateFrame}: }{\code{tkwin} class object; a special frame for some GUI parts.}
#' \item{\code{vbbox1}: }{\code{variableboxes} class object; the frame to select variables.}
#' \item{\code{vbbox2}: }{\code{variableboxes} class object; the frame to select facet variables.}
#' \item{\code{lbbox1}: }{\code{textfields} class object; the frame to set axis labels, the legend label, and the main title.}
#' \item{\code{rbbox1}: }{\code{radioboxes} class object; the frame to set the plot type.}
#' \item{\code{rbbox2}: }{\code{radioboxes} class object; the frame to set the summary type.}
#' \item{\code{rbbox3}: }{\code{radioboxes} class object; the frame to set the scale type.}
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
#' @seealso \code{\link[ggplot2:scale_date]{scale_date}}
#' @name gline-class
#' @aliases gline
#' @rdname plot-gline
#' @docType class
#' @keywords hplot
#' @importFrom tcltk2 tk2tip
#' @importFrom scales date_format
#' @export gline
gline <- setRefClass(

  Class = "gline",

  fields = c("vbbox1", "vbbox2", "lbbox1", "rbbox1", "rbbox2", "rbbox3", "tbbox1"),

  contains = c("plot_base"),

  methods = list(

    setFront = function() {

      vbbox1 <<- variableboxes$new()
      vbbox1$front(
        top       = top, 
        types     = list(nonFactors(), nonFactors(), Factors()),
        titles    = list(
          gettextKmg2("X variable (pick one)"),
          gettextKmg2("Y variable (pick one)"),
          gettextKmg2("Stratum variable")
        ),
        initialSelection = list(0, 0, FALSE)
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
        top    = alternateFrame,
        labels = list(
          gettextKmg2("Line"),
          gettextKmg2("Step"),
          gettextKmg2("Area")
        ),
        title  = gettextKmg2("Plot type")
      )

      rbbox2 <<- radioboxes$new()
      rbbox2$front(
        top    = alternateFrame,
        labels = list(
          gettextKmg2("None"),
          gettextKmg2("Mean"),
          gettextKmg2("Sum")
        ),
        title  = gettextKmg2("Summarization")
      )

      rbbox3 <<- radioboxes$new()
      rbbox3$front(
        top    = alternateFrame,
        labels = list(
          gettextKmg2("None"),
          gettextKmg2("%Y-%m-%d (ISO 8601)"),
          gettextKmg2("%m/%d/%y (ISO C99)")
        ),
        title  = gettextKmg2("Scale")
      )
      tk2tip(
        rbbox3$frame,
        gettextKmg2("If you needs complementary informations about scales, \nplease run or search ?format.Date, ?scale_x_date, or ?strptime.")
      )

      tbbox1 <<- toolbox$new()
      tbbox1$front(top)

    },

    setBack = function() {

      vbbox1$back()
      vbbox2$back()
      lbbox1$back()
      tkgrid(
        rbbox1$frame,
        labelRcmdr(alternateFrame, text = "    "),
        rbbox2$frame,
        labelRcmdr(alternateFrame, text = "    "),
        rbbox3$frame, stick = "nw")
      tkgrid(alternateFrame, stick = "nw")
      tkgrid(labelRcmdr(alternateFrame, text = "    "), stick = "nw")
      tbbox1$back()

    },

    getWindowTitle = function() {
      
      gettextKmg2("Line chart")
      
    },
    
    getHelp = function() {
      
      "geom_line"
      
    },

    getParms = function() {

      x      <- getSelection(vbbox1$variable[[1]])
      y      <- getSelection(vbbox1$variable[[2]])
      z      <- getSelection(vbbox1$variable[[3]])

      s      <- getSelection(vbbox2$variable[[1]])
      t      <- getSelection(vbbox2$variable[[2]])

      x      <- checkVariable(x)
      y      <- checkVariable(y)
      z      <- checkVariable(z)
      s      <- checkVariable(s)
      t      <- checkVariable(t)

      xlab   <- tclvalue(lbbox1$fields[[1]]$value)
      xauto  <- x
      ylab   <- tclvalue(lbbox1$fields[[2]]$value)
      yauto  <- y
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
      
      plotType    <- tclvalue(rbbox1$value)
      summaryType <- tclvalue(rbbox2$value)
      scaleType   <- tclvalue(rbbox3$value)

      list(
        x = x, y = y, z = z, s = s, t = t,
        xlab = xlab, xauto = xauto, ylab = ylab, yauto = yauto, zlab = zlab, main = main,
        size = size, family = family, colour = colour, save = save, theme = theme,
        plotType = plotType, summaryType = summaryType, scaleType = scaleType
      )

    },

    checkError = function(parms) {

      if (length(parms$x) == 0) {
        errorCondition(
          recall  = windowLine,
          message = gettextKmg2("X variable is not selected")
        )
        errorCode <- TRUE
      } else if (length(parms$y) == 0) {
        errorCondition(
          recall  = windowLine,
          message = gettextKmg2("Y variable is not selected")
        )
        errorCode <- TRUE
      } else {
        errorCode <- FALSE
      }
      errorCode

    },

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
      command <- paste0(
        ".df <- data.frame(",
        command,
        ")\n",
        ".df <- .df[order(.df$x), ]"
      )

      commandDoIt(command)
      registRmlist(.df)

    },

    getGgplot = function(parms) {

      if (length(parms$z) == 0) {
        ggplot <-  "ggplot(data = .df, aes(x = x, y = y)) + "
      } else if (any(parms$plotType == c("1", "2"))) {
        ggplot <-  "ggplot(data = .df, aes(x = x, y = y, colour = z, shape = z)) + "
      } else {
        ggplot <-  "ggplot(data = .df, aes(x = x, y = y, fill = z)) + "
      }
      ggplot

    },

    getGeom = function(parms) {

      if (parms$summaryType == "1") {
        if (parms$plotType == "1") {
          geom <- "geom_point() + geom_line(size = 1) + "
        } else if (parms$plotType == "2") {
          geom <- "geom_point() + geom_step(size = 1) + "
        } else {
          geom <- "geom_area(alpha = 0.3) + "
        }
      } else if (parms$summaryType == "2") {
        if (parms$plotType == "1") {
          geom <- "stat_summary(fun.y = \"mean\", geom = \"point\") + stat_summary(fun.y = \"mean\", geom = \"line\", size = 1) + "
        } else if (parms$plotType == "2") {
          geom <- "stat_summary(fun.y = \"mean\", geom = \"point\") + stat_summary(fun.y = \"mean\", geom = \"step\", size = 1) + "
        } else {
          geom <- "stat_summary(fun.y = \"mean\", geom = \"area\", alpha = 0.3) + "
        }
      } else {
        if (parms$plotType == "1") {
          geom <- "stat_summary(fun.y = \"sum\", geom = \"point\") + stat_summary(fun.y = \"sum\", geom = \"line\", size = 1) + "
        } else if (parms$plotType == "2") {
          geom <- "stat_summary(fun.y = \"sum\", geom = \"point\") + stat_summary(fun.y = \"sum\", geom = \"step\", size = 1) + "
        } else {
          geom <- "stat_summary(fun.y = \"sum\", geom = \"area\", alpha = 0.3) + "
        }
      }
      geom

    },

    getScale = function(parms) {
      
      if (length(parms$z) == 0) {
        scale <- ""
      } else if (any(parms$plotType == c("1", "2"))) {
        scale <- paste("scale_colour_brewer(palette = \"", parms$colour, "\") + ", sep="")
      } else {
        scale <- paste("scale_fill_brewer(palette = \"", parms$colour, "\") + ", sep="")
      }

      if (parms$scaleType == "2") {
        scale <- paste0(
          scale,
          "scale_x_date(labels = scales::date_format()) + "
        )
      } else if (parms$scaleType == "3") {
        scale <- paste0(
          scale,
          "scale_x_date(labels = scales::date_format(\"%m/%d/%y\")) + "
        )
      }

      scale <- paste0(
        scale,
        "scale_y_continuous(expand = c(0.01, 0)) + "
      )
      scale
      
    },

    getZlab = function(parms) {

      if (length(parms$z) == 0) {
        zlab <- ""
      } else if (nchar(parms$zlab) == 0) {
        zlab <- ""
      } else if (parms$zlab == "<auto>") {
        if (any(parms$plotType == c("1", "2"))) {
          zlab <- paste0("labs(colour = \"", parms$z, "\", shape = \"", parms$z, "\") + ")
        } else {
          zlab <- paste0("labs(fill = \"", parms$z, "\") + ")
        }
      } else {
        if (any(parms$plotType == c("1", "2"))) {
          zlab <- paste0("labs(colour = \"", parms$zlab, "\", shape = \"", parms$zlab, "\") + ")
        } else {
          zlab <- paste0("labs(fill = \"", parms$zlab, "\") + ")
        }
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



#' Wrapper Function of Line Chart Subclass
#'
#' \code{windowScatter} function is a wrapper function of \code{gline} class for the R-commander menu bar.
#'
#' @rdname plot-gline-windowLine
#' @keywords hplot
#' @export
windowLine <- function() {

  Line <- RcmdrPlugin.KMggplot2::gline$new()
  Line$plotWindow()

}
