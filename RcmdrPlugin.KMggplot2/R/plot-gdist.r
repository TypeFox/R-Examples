#' Distribution Plot Subclass
#'
#' \code{gdist} class is a subclass for distribution plots.
#'
#' This class is a subclass which show dialog boxes of distribution plots for graphics editing.
#'
#' @section Fields:
#' \describe{
#' \item{\code{top}: }{\code{tkwin} class object; parent of widget window.}
#' \item{\code{alternateFrame}: }{\code{tkwin} class object; a special frame for some GUI parts.}
#' \item{\code{lbbox1}: }{\code{textfields} class object; the frame to set distribution parameters.}
#' \item{\code{lbbox2}: }{\code{textfields} class object; the frame to set axis labels, the legend label, and the main title.}
#' \item{\code{rbbox1}: }{\code{radioboxes} class object; the frame to set the function type.}
#' \item{\code{tbbox1}: }{\code{toolbox} class object; the frame to set the font, the colour set, other option, and the theme.}
#' \item{\code{windowTitle}: }{Character; the window title.}
#' \item{\code{distType}: }{Character; the distribution type ("discrete" or "continuous").}
#' \item{\code{distName}: }{Character; the distribution name.}
#' \item{\code{parmNames}: }{List of Characters; names of distribution parameters.}
#' \item{\code{parmInits}: }{List of Characters; initial values of distribution parameters.}
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
#' @name gdist-class
#' @aliases gdist
#' @rdname plot-gdist
#' @docType class
#' @keywords hplot
#' @export gdist
gdist <- setRefClass(

  Class = "gdist",

  fields = c("lbbox1", "lbbox2", "rbbox1", "tbbox1", "windowTitle", "distType", "distName", "parmNames", "parmInits"),

  contains = c("plot_base"),

  methods = list(

    setFront = function() {

      lbbox1 <<- textfields$new()
      lbbox1$front(
        top        = top,
        initValues = parmInits,
        titles     = parmNames
      )

      lbbox2 <<- textfields$new()
      lbbox2$front(
        top        = top,
        initValues = list("<auto>", "<auto>", "<auto>"),
        titles     = list(
          gettextKmg2("Horizontal axis label"),
          gettextKmg2("Vertical axis label"),
          gettextKmg2("Title")
        )
      )

      rbbox1 <<- radioboxes$new()
      rbbox1$front(
        top    = top,
        labels = list(
          gettextKmg2("Plot density function"),
          gettextKmg2("Plot distribution function")
        ),
        title  = gettextKmg2("Function type")
      )

      tbbox1 <<- toolbox$new()
      tbbox1$front(top, showcolourbox = FALSE)

    },

    setBack = function() {

      lbbox1$back()
      lbbox2$back()
      rbbox1$back()
      tbbox1$back(4)

    },

    getWindowTitle = function() {
      
      windowTitle
      
    },
    
    getHelp = function() {
      
      "Distributions"
      
    },

    getParms = function() {

      x      <- character(0)
      y      <- character(0)
      z      <- character(0)

      s      <- character(0)
      t      <- character(0)

      xlab   <- tclvalue(lbbox2$fields[[1]]$value)
      xauto  <- "x"
      ylab   <- tclvalue(lbbox2$fields[[2]]$value)
      # yauto  <- y
      zlab   <- character(0)
      main   <- tclvalue(lbbox2$fields[[3]]$value)

      size   <- tclvalue(tbbox1$size$value)
      family <- getSelection(tbbox1$family)
      colour <- character(0)
      save   <- tclvalue(tbbox1$goption$value[[1]])
      theme  <- checkTheme(getSelection(tbbox1$theme))
      
      options(
        kmg2FontSize   = tclvalue(tbbox1$size$value),
        kmg2FontFamily = seq_along(tbbox1$family$varlist)[tbbox1$family$varlist == getSelection(tbbox1$family)] - 1,
        kmg2SaveGraph  = tclvalue(tbbox1$goption$value[[1]]),
        kmg2Theme      = seq_along(tbbox1$theme$varlist)[tbbox1$theme$varlist == getSelection(tbbox1$theme)] - 1
      )
      
      funcType <- tclvalue(rbbox1$value)
      if (funcType == "1") {
        yauto <- "Density"
      } else {
        yauto <- "Cumulative Probability"
      }
      
      parmLength  <- length(parmInits)
      parmValues  <- lapply(1:parmLength,
        function(i, lbbox1) tclvalue(lbbox1$fields[[i]]$value), lbbox1)
      parmValuesList <- ""
      for (i in 1:parmLength) {
        parmValuesList <- paste0(parmValuesList, ", ", parmValues[i])
      }

      list(
        x = x, y = y, z = z, s = s, t = t,
        xlab = xlab, xauto = xauto, ylab = ylab, yauto = yauto, zlab = zlab, main = main,
        size = size, family = family, colour = colour, save = save, theme = theme,
        funcType = funcType, parmValuesList = parmValuesList
      )

    },

    setDataframe = function(parms) {

      command <- paste0("q", distName, "(c(0.001, 1-0.001)", parms$parmValuesList, ")")
      range <- eval(parse(text = command))

      if (distType == "continuous") {
        command <- paste0(".x  <- seq(", range[1], ", ", range[2], ", length.out = 100)")
      } else {
        command <- paste0(".x  <- ", range[1], ":", range[2])
      }
      commandDoIt(command)
      registRmlist(.x)
      
      if (parms$funcType == "1") {
        command <- paste0(".df <- data.frame(x = .x, y = d", distName, "(.x", parms$parmValuesList, "))")
      } else {
        command <- paste0(".df <- data.frame(x = .x, y = p", distName, "(.x", parms$parmValuesList, "))")
      }
      commandDoIt(command)
      registRmlist(.df)

    },

    getGgplot = function(parms) {

      "ggplot(.df, aes(x = x, y = y)) + "

    },

    getGeom = function(parms) {
      
      if (distType == "continuous") {
        geom <- "geom_line(size = 1.5) + "
      } else if (distType == "discrete" && parms$funcType == "1") {
        geom <- "geom_bar(stat = \"identity\") + "
      } else {
        geom <- "geom_step(size = 1.5) + "
      }
      geom

    },

    getMain = function(parms) {

      if (nchar(parms$main) == 0) {
        main <- ""
      } else if (parms$main == "<auto>") {
        if (parms$funcType == "1") {
          main <- paste0("labs(title = \"d", distName, "(x", parms$parmValuesList, ")\") + ")
        } else {
          main <- paste0("labs(title = \"p", distName, "(x", parms$parmValuesList, ")\") + ")
        }
      } else {
        main <- paste0("labs(title = \"", parms$main, "\") + ")
      }
      main

    }

  )
)



#' Wrapper Function of Normal Distribution Plot Subclass
#'
#' \code{windowDistNorm} function is a wrapper function of \code{gdist} class for the R-commander menu bar.
#'
#' @rdname plot-gdist-windowDistNorm
#' @keywords hplot
#' @export
windowDistNorm <- function() {

  DistNorm <- RcmdrPlugin.KMggplot2::gdist$new(
    windowTitle = gettextKmg2("Plot normal distiribution"),
    distType    = "continuous",
    distName    = "norm",
    parmNames  = list(
      gettextKmg2("Mean"),
      gettextKmg2("S.D.")
    ),
    parmInits   = list("0", "1")
  )
  DistNorm$plotWindow()

}



#' Wrapper Function of t Distribution Plot Subclass
#'
#' \code{windowDistT} function is a wrapper function of \code{gdist} class for the R-commander menu bar.
#'
#' @rdname plot-gdist-windowDistT
#' @keywords hplot
#' @export
windowDistT <- function() {

  DistT <- RcmdrPlugin.KMggplot2::gdist$new(
    windowTitle = gettextKmg2("Plot t distiribution"),
    distType    = "continuous",
    distName    = "t",
    parmNames  = list(
      gettextKmg2("df")
    ),
    parmInits   = list("5")
  )
  DistT$plotWindow()

}



#' Wrapper Function of Chi-square Distribution Plot Subclass
#'
#' \code{windowDistChisq} function is a wrapper function of \code{gdist} class for the R-commander menu bar.
#'
#' @rdname plot-gdist-windowDistChisq
#' @keywords hplot
#' @export
windowDistChisq <- function() {

  DistChisq <- RcmdrPlugin.KMggplot2::gdist$new(
    windowTitle = gettextKmg2("Plot chi-square distiribution"),
    distType    = "continuous",
    distName    = "chisq",
    parmNames  = list(
      gettextKmg2("df")
    ),
    parmInits   = list("5")
  )
  DistChisq$plotWindow()

}



#' Wrapper Function of F Distribution Plot Subclass
#'
#' \code{windowDistF} function is a wrapper function of \code{gdist} class for the R-commander menu bar.
#'
#' @rdname plot-gdist-windowDistF
#' @keywords hplot
#' @export
windowDistF <- function() {

  DistF <- RcmdrPlugin.KMggplot2::gdist$new(
    windowTitle = gettextKmg2("Plot F distiribution"),
    distType    = "continuous",
    distName    = "f",
    parmNames  = list(
      gettextKmg2("Numerator df"),
      gettextKmg2("Denominator df")
    ),
    parmInits   = list("2", "3")
  )
  DistF$plotWindow()

}



#' Wrapper Function of Exponential Distribution Plot Subclass
#'
#' \code{windowDistExp} function is a wrapper function of \code{gdist} class for the R-commander menu bar.
#'
#' @rdname plot-gdist-windowDistExp
#' @keywords hplot
#' @export
windowDistExp <- function() {

  DistExp <- RcmdrPlugin.KMggplot2::gdist$new(
    windowTitle = gettextKmg2("Plot exponential distiribution"),
    distType    = "continuous",
    distName    = "exp",
    parmNames  = list(
      gettextKmg2("rate")
    ),
    parmInits   = list("1")
  )
  DistExp$plotWindow()

}



#' Wrapper Function of Uniform Distribution Plot Subclass
#'
#' \code{windowDistUnif} function is a wrapper function of \code{gdist} class for the R-commander menu bar.
#'
#' @rdname plot-gdist-windowDistUnif
#' @keywords hplot
#' @export
windowDistUnif <- function() {

  DistUnif <- RcmdrPlugin.KMggplot2::gdist$new(
    windowTitle = gettextKmg2("Plot uniform distiribution"),
    distType    = "continuous",
    distName    = "unif",
    parmNames  = list(
      gettextKmg2("Minimum"),
      gettextKmg2("Maximum")
    ),
    parmInits   = list("0", "1")
  )
  DistUnif$plotWindow()

}



#' Wrapper Function of Beta Distribution Plot Subclass
#'
#' \code{windowDistBeta} function is a wrapper function of \code{gdist} class for the R-commander menu bar.
#'
#' @rdname plot-gdist-windowDistBeta
#' @keywords hplot
#' @export
windowDistBeta <- function() {

  DistBeta <- RcmdrPlugin.KMggplot2::gdist$new(
    windowTitle = gettextKmg2("Plot beta distiribution"),
    distType    = "continuous",
    distName    = "beta",
    parmNames  = list(
      gettextKmg2("Shape 1"),
      gettextKmg2("Shape 2")
    ),
    parmInits   = list("9", "3")
  )
  DistBeta$plotWindow()

}



#' Wrapper Function of Cauchy Distribution Plot Subclass
#'
#' \code{windowDistCauchy} function is a wrapper function of \code{gdist} class for the R-commander menu bar.
#'
#' @rdname plot-gdist-windowDistCauchy
#' @keywords hplot
#' @export
windowDistCauchy <- function() {

  DistCauchy <- RcmdrPlugin.KMggplot2::gdist$new(
    windowTitle = gettextKmg2("Plot cauchy distiribution"),
    distType    = "continuous",
    distName    = "cauchy",
    parmNames  = list(
      gettextKmg2("Location"),
      gettextKmg2("Scale")
    ),
    parmInits   = list("0", "1")
  )
  DistCauchy$plotWindow()

}



#' Wrapper Function of Logistic Distribution Plot Subclass
#'
#' \code{windowDistLogis} function is a wrapper function of \code{gdist} class for the R-commander menu bar.
#'
#' @rdname plot-gdist-windowDistLogis
#' @keywords hplot
#' @export
windowDistLogis <- function() {

  DistLogis <- RcmdrPlugin.KMggplot2::gdist$new(
    windowTitle = gettextKmg2("Plot logistic distiribution"),
    distType    = "continuous",
    distName    = "logis",
    parmNames  = list(
      gettextKmg2("Location"),
      gettextKmg2("Scale")
    ),
    parmInits   = list("0", "1")
  )
  DistLogis$plotWindow()

}



#' Wrapper Function of Log-normal Distribution Plot Subclass
#'
#' \code{windowDistLnorm} function is a wrapper function of \code{gdist} class for the R-commander menu bar.
#'
#' @rdname plot-gdist-windowDistLnorm
#' @keywords hplot
#' @export
windowDistLnorm <- function() {

  DistLnorm <- RcmdrPlugin.KMggplot2::gdist$new(
    windowTitle = gettextKmg2("Plot log-normal distiribution"),
    distType    = "continuous",
    distName    = "lnorm",
    parmNames  = list(
      gettextKmg2("Mean (log scale)"),
      gettextKmg2("S.D. (log scale)")
    ),
    parmInits   = list("0", "1")
  )
  DistLnorm$plotWindow()

}



#' Wrapper Function of Gamma Distribution Plot Subclass
#'
#' \code{windowDistGamma} function is a wrapper function of \code{gdist} class for the R-commander menu bar.
#'
#' @rdname plot-gdist-windowDistGamma
#' @keywords hplot
#' @export
windowDistGamma <- function() {

  DistGamma <- RcmdrPlugin.KMggplot2::gdist$new(
    windowTitle = gettextKmg2("Plot gamma distiribution"),
    distType    = "continuous",
    distName    = "gamma",
    parmNames  = list(
      gettextKmg2("Shape"),
      gettextKmg2("Rate (inverse scale)")
    ),
    parmInits   = list("1", "1")
  )
  DistGamma$plotWindow()

}



#' Wrapper Function of Weibull Distribution Plot Subclass
#'
#' \code{windowDistWeibull} function is a wrapper function of \code{gdist} class for the R-commander menu bar.
#'
#' @rdname plot-gdist-windowDistWeibull
#' @keywords hplot
#' @export
windowDistWeibull <- function() {

  DistWeibull <- RcmdrPlugin.KMggplot2::gdist$new(
    windowTitle = gettextKmg2("Plot weibull distiribution"),
    distType    = "continuous",
    distName    = "weibull",
    parmNames  = list(
      gettextKmg2("Shape"),
      gettextKmg2("Scale")
    ),
    parmInits   = list("1", "pi")
  )
  DistWeibull$plotWindow()

}



#' Wrapper Function of Binomial Distribution Plot Subclass
#'
#' \code{windowDistBinom} function is a wrapper function of \code{gdist} class for the R-commander menu bar.
#'
#' @rdname plot-gdist-windowDistBinom
#' @keywords hplot
#' @export
windowDistBinom <- function() {

  DistBinom <- RcmdrPlugin.KMggplot2::gdist$new(
    windowTitle = gettextKmg2("Plot binomial distiribution"),
    distType    = "discrete",
    distName    = "binom",
    parmNames  = list(
      gettextKmg2("Binomial trials"),
      gettextKmg2("Probability of success")
    ),
    parmInits   = list("20", "0.5")
  )
  DistBinom$plotWindow()

}



#' Wrapper Function of Poisson Distribution Plot Subclass
#'
#' \code{windowDistPois} function is a wrapper function of \code{gdist} class for the R-commander menu bar.
#'
#' @rdname plot-gdist-windowDistPois
#' @keywords hplot
#' @export
windowDistPois <- function() {

  DistPois <- RcmdrPlugin.KMggplot2::gdist$new(
    windowTitle = gettextKmg2("Plot poisson distiribution"),
    distType    = "discrete",
    distName    = "pois",
    parmNames  = list(
      gettextKmg2("Mean")
    ),
    parmInits   = list("10")
  )
  DistPois$plotWindow()

}



#' Wrapper Function of Geometric Distribution Plot Subclass
#'
#' \code{windowDistGeom} function is a wrapper function of \code{gdist} class for the R-commander menu bar.
#'
#' @rdname plot-gdist-windowDistGeom
#' @keywords hplot
#' @export
windowDistGeom <- function() {

  DistGeom <- RcmdrPlugin.KMggplot2::gdist$new(
    windowTitle = gettextKmg2("Plot geometric distiribution"),
    distType    = "discrete",
    distName    = "geom",
    parmNames  = list(
      gettextKmg2("Probability of success")
    ),
    parmInits   = list("0.25")
  )
  DistGeom$plotWindow()

}



#' Wrapper Function of Hypergeometric Distribution Plot Subclass
#'
#' \code{windowDistHyper} function is a wrapper function of \code{gdist} class for the R-commander menu bar.
#'
#' @rdname plot-gdist-windowDistHyper
#' @keywords hplot
#' @export
windowDistHyper <- function() {

  DistHyper <- RcmdrPlugin.KMggplot2::gdist$new(
    windowTitle = gettextKmg2("Plot hypergeometric distiribution"),
    distType    = "discrete",
    distName    = "hyper",
    parmNames  = list(
      gettextKmg2("m (number of white balls in the urn)"),
      gettextKmg2("n (number of black balls in the urn)"),
      gettextKmg2("k (number of balls drawn from the urn)")
    ),
    parmInits   = list("25", "20", "8")
  )
  DistHyper$plotWindow()

}



#' Wrapper Function of Negative Binomial Distribution Plot Subclass
#'
#' \code{windowDistNbinom} function is a wrapper function of \code{gdist} class for the R-commander menu bar.
#'
#' @rdname plot-gdist-windowDistNbinom
#' @keywords hplot
#' @export
windowDistNbinom <- function() {

  DistNbinom <- RcmdrPlugin.KMggplot2::gdist$new(
    windowTitle = gettextKmg2("Plot negative binomial distiribution"),
    distType    = "discrete",
    distName    = "nbinom",
    parmNames  = list(
      gettextKmg2("Target number of success"),
      gettextKmg2("Probability of success")
    ),
    parmInits   = list("5", "0.5")
  )
  DistNbinom$plotWindow()

}
