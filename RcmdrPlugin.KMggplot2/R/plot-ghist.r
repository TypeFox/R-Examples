#' Histogram Subclass
#'
#' \code{ghist} class is a subclass for histograms.
#'
#' This class is a subclass which show dialog boxes of histograms for graphics editing.
#'
#' @section Fields:
#' \describe{
#' \item{\code{top}: }{\code{tkwin} class object; parent of widget window.} 
#' \item{\code{alternateFrame}: }{\code{tkwin} class object; a special frame for some GUI parts.} 
#' \item{\code{vbbox1}: }{\code{variableboxes} class object; the frame to select variables.} 
#' \item{\code{vbbox2}: }{\code{variableboxes} class object; the frame to select facet variables.} 
#' \item{\code{vbbox3}: }{\code{variableboxes} class object; the frame to set the no. of bins.} 
#' \item{\code{lbbox1}: }{\code{textfields} class object; the frame to set axis labels and the main title.} 
#' \item{\code{rbbox1}: }{\code{radioboxes} class object; the frame to set the axis scaling.} 
#' \item{\code{cbbox1}: }{\code{checkboxes} class object; the frame to set options.} 
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
#' @name ghist-class
#' @aliases ghist
#' @rdname plot-ghist
#' @docType class
#' @keywords hplot
#' @importFrom scales percent_format
#' @importFrom RColorBrewer brewer.pal
#' @export ghist
ghist <- setRefClass(

  Class = "ghist",

  fields = c("vbbox1", "vbbox2", "vbbox3", "lbbox1", "rbbox1", "cbbox1", "tbbox1"),

  contains = c("plot_base"),

  methods = list(

    setFront = function() {

      vbbox1 <<- variableboxes$new()
      vbbox1$front(
        top       = top, 
        types     = list(nonFactors()),
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
        initValues = list("<auto>", "<auto>", ""),
        titles     = list(
          gettextKmg2("Horizontal axis label"),
          gettextKmg2("Vertical axis label"),
          gettextKmg2("Title")
        )
      )

      vbbox3 <<- variableboxes$new()
      vbbox3$front(
        top       = alternateFrame, 
        types     = list(c("Scott", "Freedman-Diaconis", "Sturges")),
        titles    = list(
          gettextKmg2("No. of bins")
        ),
        initialSelection = list(0)
      )

      rbbox1 <<- radioboxes$new()
      rbbox1$front(
        top    = alternateFrame,
        labels = list(
          gettextKmg2("Densities"),
          gettextKmg2("Frequency counts"),
          gettextKmg2("Percentages")
        ),
        title  = gettextKmg2("Axis scaling")
      )

      cbbox1 <<- checkboxes$new()
      cbbox1$front(
        top        = alternateFrame,
        initValues = list("0", "0"),
        labels     = list(
          gettextKmg2("Density estimation"),
          gettextKmg2("Heat map")
        ),
        title      = gettextKmg2("Options")
      )

      tbbox1 <<- toolbox$new()
      tbbox1$front(top)

    },

    setBack = function() {

      vbbox1$back()
      vbbox2$back()
      lbbox1$back()

      boxlist <- c(
        list(vbbox3$frame),
        list(labelRcmdr(alternateFrame, text="    ")),
        list(cbbox1$frame),
        list(labelRcmdr(alternateFrame, text="    ")),
        list(rbbox1$frame)
      )
      do.call(tkgrid, c(vbbox3$back_list, list(sticky="nw")))
      do.call(tkgrid, c(boxlist, list(sticky="nw")))
      tkgrid(alternateFrame, stick="nw")
      tkgrid(labelRcmdr(alternateFrame, text="    "), stick="nw")

      tbbox1$back()

    },

    getWindowTitle = function() {
      
      gettextKmg2("Histogram")
      
    },
    
    getHelp = function() {
      
      "geom_histogram"
      
    },

    getParms = function() {

      x      <- getSelection(vbbox1$variable[[1]])
      # y      <- ""
      z      <- character(0)

      s      <- getSelection(vbbox2$variable[[1]])
      t      <- getSelection(vbbox2$variable[[2]])

      x      <- checkVariable(x)
      s      <- checkVariable(s)
      t      <- checkVariable(t)

      xlab   <- tclvalue(lbbox1$fields[[1]]$value)
      xauto  <- x
      ylab   <- tclvalue(lbbox1$fields[[2]]$value)
      # yauto  <- y
      zlab   <- ""
      main   <- tclvalue(lbbox1$fields[[3]]$value)

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
      
      densityPlot  <- tclvalue(cbbox1$value[[1]])
      heatPlot     <- tclvalue(cbbox1$value[[2]])
      nbins        <- getSelection(vbbox3$variable[[1]])
      axisScaling  <- tclvalue(rbbox1$value)
      if (densityPlot == "1" || axisScaling == "1") {
        yauto       <- gettextKmg2("Density")
        y           <- "..density.."
        axisScaling <- "1"
      }  else if (axisScaling == "2") {
        yauto  <- gettextKmg2("Count")
        y      <- "..count.."
      } else if (axisScaling == "3") {
        yauto  <- gettextKmg2("Percent")
        y      <- "..count../sum(..count..)"
      }

      list(
        x = x, y = y, z = z, s = s, t = t,
        xlab = xlab, xauto = xauto, ylab = ylab, yauto = yauto, zlab = zlab, main = main,
        size = size, family = family, colour = colour, save = save, theme = theme,
        axisScaling = axisScaling, densityPlot = densityPlot, heatPlot = heatPlot, nbins = nbins
      )

    },

    checkError = function(parms) {

      if (length(parms$x) == 0) {
        errorCondition(
          recall  = windowHist,
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
      if (length(parms$x) != 0) {
        var <- c(var, paste0("x = ", ActiveDataSet(), "$", parms$x))
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

    getGgplot = function(parms) {

      paste0(
        "ggplot(data = .df, aes(x = x, y = ", parms$y, ")) + "
      )

    },

    getGeom = function(parms) {

      if (length(parms$nbins) == 0) {
        command <- ".nbins <- pretty(range(.df$x), n = nclass.scott(.df$x), min.n = 1)"
      } else if (parms$nbins == "Sturges") {
        command <- ".nbins <- pretty(range(.df$x), n = nclass.Sturges(.df$x), min.n = 1)"
      } else if (parms$nbins == "Freedman-Diaconis") {
        command <- ".nbins <- pretty(range(.df$x), n = nclass.FD(.df$x), min.n = 1)"
      } else {
        command <- ".nbins <- pretty(range(.df$x), n = nclass.scott(.df$x), min.n = 1)"
      }
      commandDoIt(command)
      registRmlist(.nbins)

      if (parms$heatPlot == "1") {
        geom <- paste0(
          "geom_histogram(aes(fill = ", parms$y, "), breaks = .nbins) + "
        )
      } else {
        geom <- "geom_histogram(breaks = .nbins) + "
      }

      if (parms$densityPlot == "1") {
        geom <- paste0(
          geom,
          "stat_density(geom = \"path\", size = 1, alpha = 0.5) + "
        )
      }
      geom

    },

    getScale = function(parms) {
      
      if (parms$axisScaling == "3") {
        scale <- "scale_y_continuous(expand = c(0.01, 0), labels = scales::percent_format()) + "
      } else {
        scale <- "scale_y_continuous(expand = c(0.01, 0)) + "
      }

      if (parms$heatPlot == "1") {
        if (parms$axisScaling == "3") {
          scale <- paste0(
            scale,
            "scale_fill_gradient(",
              "low = RColorBrewer::brewer.pal(3, \"", parms$colour, "\")[2], ", 
              "high = RColorBrewer::brewer.pal(3, \"", parms$colour, "\")[1], ",
              "labels = scales::percent_format()",
            ") + "
          )
        } else {
          scale <- paste0(
            scale,
            "scale_fill_gradient(",
              "low = RColorBrewer::brewer.pal(3, \"", parms$colour, "\")[2], ", 
              "high = RColorBrewer::brewer.pal(3, \"", parms$colour, "\")[1]",
            ") + "
          )
        }
      }
      scale

    },

    getOpts = function(parms) {

      opts <- list()
      if (length(parms$s) != 0 || length(parms$t) != 0) {
        opts <- c(opts, "panel.margin = unit(0.3, \"lines\")")
      }

      if (parms$heatPlot == "1") {
        opts <- c(opts, "legend.position = \"none\"")
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



#' Wrapper Function of Histogram Subclass
#'
#' \code{windowHist} function is a wrapper function of \code{ghist} class for the R-commander menu bar.
#'
#' @rdname plot-ghist-windowHist
#' @keywords hplot
#' @export
windowHist <- function() {

  Hist <- RcmdrPlugin.KMggplot2::ghist$new()
  Hist$plotWindow()

}
