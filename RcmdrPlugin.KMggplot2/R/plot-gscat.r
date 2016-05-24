#' Scatter Plot Subclass
#'
#' \code{gscat} class is a subclass for scatter plots.
#'
#' This class is a subclass which show dialog boxes of scatter plots for graphics editing.
#'
#' @section Fields:
#' \describe{
#' \item{\code{top}: }{\code{tkwin} class object; parent of widget window.}
#' \item{\code{alternateFrame}: }{\code{tkwin} class object; a special frame for some GUI parts.}
#' \item{\code{vbbox1}: }{\code{variableboxes} class object; the frame to select variables.}
#' \item{\code{vbbox2}: }{\code{variableboxes} class object; the frame to select facet variables.}
#' \item{\code{lbbox1}: }{\code{textfields} class object; the frame to set axis labels, the legend label, and the main title.}
#' \item{\code{rbbox1}: }{\code{radioboxes} class object; the frame to set the smoothing type.}
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
#' @name gscat-class
#' @aliases gscat
#' @rdname plot-gscat
#' @docType class
#' @keywords hplot
#' @export gscat
gscat <- setRefClass(

  Class = "gscat",

  fields = c("vbbox1", "vbbox2", "lbbox1", "rbbox1", "tbbox1"),

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
        top    = top,
        labels = list(
          gettextKmg2("None"),
          gettextKmg2("Smoothing with C.I. (linear regression)"),
          gettextKmg2("Smoothing without C.I. (linear regression)"),
          gettextKmg2("Smoothing with C.I. (loess or gam)"),
          gettextKmg2("Smoothing without C.I. (loess or gam)")
        ),
        title  = gettextKmg2("Smoothing type")
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
      
      gettextKmg2("Scatter plot")
      
    },
    
    getHelp = function() {
      
      "stat_smooth"
      
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
      
      smoothType   <- tclvalue(rbbox1$value)

      list(
        x = x, y = y, z = z, s = s, t = t,
        xlab = xlab, xauto = xauto, ylab = ylab, yauto = yauto, zlab = zlab, main = main,
        size = size, family = family, colour = colour, save = save, theme = theme,
        smoothType = smoothType
      )

    },

    checkError = function(parms) {

      if (length(parms$x) == 0) {
        errorCondition(
          recall  = windowScatter,
          message = gettextKmg2("X variable is not selected")
        )
        errorCode <- TRUE
      } else if (length(parms$y) == 0) {
        errorCondition(
          recall  = windowScatter,
          message = gettextKmg2("Y variable is not selected")
        )
        errorCode <- TRUE
      } else {
        errorCode <- FALSE
      }
      errorCode

    },

    getGgplot = function(parms) {

      if (length(parms$z) == 0) {
        ggplot <-  "ggplot(data = .df, aes(x = x, y = y)) + "
      } else {
        ggplot <-  "ggplot(data = .df, aes(x = x, y = y, colour = z, shape = z)) + "
      }
      ggplot

    },

    getGeom = function(parms) {

      geom <- "geom_point() + "

      if (length(parms$z) == 0) {
        aes <- ""
      } else if (parms$smoothType == "4") {
        aes <- "aes(fill = z)"
      } else {
        aes <- "aes(fill = z), "
      }

      if (parms$smoothType == "2") {
        geom <- paste0(
          geom,
          "stat_smooth(", aes, "method = \"lm\") + "
        )
      } else if (parms$smoothType == "3") {
        geom <- paste0(
          geom,
          "stat_smooth(", aes, "method = \"lm\", se = FALSE) + "
        )
      } else if (parms$smoothType == "4") {
        geom <- paste0(
          geom,
          "stat_smooth(", aes, ") + "
        )
      } else if (parms$smoothType == "5") {
        geom <- paste0(
          geom,
          "stat_smooth(", aes, "se = FALSE) + "
        )
      }
      geom

    },

    getScale = function(parms) {
      
      scale <- "scale_y_continuous(expand = c(0.01, 0)) + "
      if (length(parms$z) != 0) {
        scale <- paste0(
          scale,
          "scale_colour_brewer(palette = \"", parms$colour, "\") + "
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
        if (parms$smoothType == "1") {
          zlab <- paste0("labs(colour = \"", parms$z, "\", shape = \"", parms$z, "\") + ")
        } else {
          zlab <- paste0("labs(colour = \"", parms$z, "\", shape = \"", parms$z, "\", fill = \"", parms$z, "\") + ")
        }
      } else {
        if (parms$smoothType == "1") {
          zlab <- paste0("labs(colour = \"", parms$zlab, "\", shape = \"", parms$zlab, "\") + ")
        } else {
          zlab <- paste0("labs(colour = \"", parms$zlab, "\", shape = \"", parms$zlab, "\", fill = \"", parms$zlab, "\") + ")
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

    },

    getMessage = function() {

      gettextKmg2("Smoothing failed.  Please try another smoothing type, or check the data and variables.")

    }

  )
)



#' Wrapper Function of Scatter Plot Subclass
#'
#' \code{windowScatter} function is a wrapper function of \code{gscat} class for the R-commander menu bar.
#'
#' @rdname plot-gscat-windowScatter
#' @keywords hplot
#' @export
windowScatter <- function() {

  Scatter <- RcmdrPlugin.KMggplot2::gscat$new()
  Scatter$plotWindow()

}
