#' Q-Q Plot Subclass
#'
#' \code{gqq} class is a subclass for Q-Q plots.
#'
#' This class is a subclass which show dialog boxes of Q-Q plots for graphics editing.
#'
#' @section Fields:
#' \describe{
#' \item{\code{top}: }{\code{tkwin} class object; parent of widget window.}
#' \item{\code{alternateFrame}: }{\code{tkwin} class object; a special frame for some GUI parts.}
#' \item{\code{vbbox1}: }{\code{variableboxes} class object; the frame to select variables.}
#' \item{\code{lbbox1}: }{\code{textfields} class object; the frame to set axis labels, the legend label, and the main title.}
#' \item{\code{lbbox2}: }{\code{textfields} class object; the frame to set other theoretical distribution.}
#' \item{\code{rbbox1}: }{\code{radioboxes} class object; the frame to set the theoretical distribution.}
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
#' @name gqq-class
#' @aliases gqq
#' @rdname plot-gqq
#' @docType class
#' @keywords hplot
#' @export gqq
gqq <- setRefClass(

  Class = "gqq",

  fields = c("vbbox1", "vbbox2", "lbbox1", "lbbox2", "rbbox1", "tbbox1"),

  contains = c("plot_base"),

  methods = list(

    setFront = function() {

      vbbox1 <<- variableboxes$new()
      vbbox1$front(
        top       = top, 
        types     = list(nonFactors()),
        titles    = list(
          gettextKmg2("Y variable (pick one)")
        ),
        initialSelection = list(0)
      )

      lbbox1 <<- textfields$new()
      lbbox1$front(
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
        top    = alternateFrame,
        labels = list(
          gettextKmg2("Normal distribution"),
          gettextKmg2("Log-normal distribution"),
          gettextKmg2("Beta distribution"),
          gettextKmg2("Exponential distribution"),
          gettextKmg2("Gamma distribution"),
          gettextKmg2("Weibull distribution"),
          gettextKmg2("Other distribution")
        ),
        title = gettextKmg2("Distribution")
      )

      lbbox2 <<- textfields$new()
      lbbox2$front(
        top        = alternateFrame,
        initValues = list(
          "distribution = qnorm, dparams = list(mean = 0, sd = 1)"
        ),
        titles     = list(
          gettextKmg2("Parameters for other distribution")
        )
      )

      tbbox1 <<- toolbox$new()
      tbbox1$front(top, showcolourbox = FALSE)

    },

    setBack = function() {

      vbbox1$back()
      lbbox1$back()

      boxlist <- c(
        list(rbbox1$frame),
        list(labelRcmdr(alternateFrame, text="    ")),
        list(lbbox2$frame)
      )
      do.call(tkgrid, c(lbbox2$back_list, list(sticky="nw")))
      do.call(tkgrid, c(boxlist, list(sticky="nw")))
      tkgrid(alternateFrame, stick="nw")
      tkgrid(labelRcmdr(alternateFrame, text="    "), stick="nw")

      tbbox1$back(4)

    },

    getWindowTitle = function() {
      
      gettextKmg2("Q-Q plot")
      
    },
    
    getHelp = function() {
      
      "Distributions"
      
    },

    getParms = function() {

      x      <- character(0)
      y      <- getSelection(vbbox1$variable[[1]])
      z      <- character(0)

      s      <- character(0)
      t      <- character(0)

      y      <- checkVariable(y)

      xlab   <- tclvalue(lbbox1$fields[[1]]$value)
      xauto  <- "Theoretical quantile"
      ylab   <- tclvalue(lbbox1$fields[[2]]$value)
      yauto  <- "Sample quantile"
      zlab   <- character(0)
      main   <- tclvalue(lbbox1$fields[[3]]$value)

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
      
      distType  <- tclvalue(rbbox1$value)
      distParms <- tclvalue(lbbox2$fields[[1]]$value)

      list(
        x = x, y = y, z = z, s = s, t = t,
        xlab = xlab, xauto = xauto, ylab = ylab, yauto = yauto, zlab = zlab, main = main,
        size = size, family = family, colour = colour, save = save, theme = theme,
        distType = distType, distParms = distParms
      )

    },

    checkError = function(parms) {

      if (length(parms$y) == 0) {
        errorCondition(
          recall  = windowQQ,
          message = gettextKmg2("Y variable is not selected")
        )
        errorCode <- TRUE
      } else {

        if (mode == 1) {
          logger("require(\"ggplot2\")")
        }

        setDataframe(parms)

        if (parms$distType == "2" && any(.df$y <= 0)) {
          response <- tclvalue(RcmdrTkmessageBox(
            message = gettextKmg2("The log-normal distribution defined on the interval (0, +Inf)."),
            title   = gettextKmg2("Error"),
            icon    = "error",
            type    = "ok",
            default = "ok")
          )
          if (response == "ok") {
            return(TRUE)
          }
        } else if (parms$distType == "3" && any(.df$y > 1 | .df$y < 0)) {
          response <- tclvalue(RcmdrTkmessageBox(
            message = gettextKmg2("The beta distribution defined on the interval [0, 1]."),
            title   = gettextKmg2("Error"),
            icon    = "error",
            type    = "ok",
            default = "ok")
          )
          if (response == "ok") {
            return(TRUE)
          }
        } else if (parms$distType == "4" && any(.df$y < 0)) {
          response <- tclvalue(RcmdrTkmessageBox(
            message = gettextKmg2("The exponential distribution defined on the interval [0, +Inf)."),
            title   = gettextKmg2("Error"),
            icon    = "error",
            type    = "ok",
            default = "ok")
          )
          if (response == "ok") {
            return(TRUE)
          }
        } else if (parms$distType == "5" && any(.df$y < 0)) {
          response <- tclvalue(RcmdrTkmessageBox(
            message = gettextKmg2("The gamma distribution defined on the interval [0, +Inf)."),
            title   = gettextKmg2("Error"),
            icon    = "error",
            type    = "ok",
            default = "ok")
          )
          if (response == "ok") {
            return(TRUE)
          }
        } else if (parms$distType == "6" && any(.df$y < 0)) {
          response <- tclvalue(RcmdrTkmessageBox(
            message = gettextKmg2("The weibull distribution defined on the interval [0, +Inf)."),
            title   = gettextKmg2("Error"),
            icon    = "error",
            type    = "ok",
            default = "ok")
          )
          if (response == "ok") {
            return(TRUE)
          }
        } else {

          if (parms$distType == "1") {
            command <- paste0(
              "# unbiased estimator\n",
              ".est <- c(mean(.df$y), sd(.df$y))"
            )
          } else if (parms$distType == "2") {
            command <- paste0(
              "# unbiased estimator\n",
              ".est <- c(mean(log(.df$y)), sd(log(.df$y)))"
            )
          } else if (parms$distType == "3") {
            command <- paste0(
              "# method-of-moments estimator\n",
              ".est <- c(mean(.df$y)*(mean(.df$y)*(1 - mean(.df$y))/var(.df$y) - 1), (1 - mean(.df$y))*(mean(.df$y)*(1 - mean(.df$y))/var(.df$y) - 1))"
            )
          } else if (parms$distType == "4") {
            command <- paste0(
              "# maximum likelihood estimator\n",
              ".est <- 1 / mean(.df$y)"
            )
          } else if (parms$distType == "5") {
            command <- paste0(
              "# maximum likelihood estimator\n",
              ".objf <- function(p, y) {",
              "-sum(dgamma(y, shape = p[1], scale = p[2], log = TRUE))",
              "}\n",
              ".est <- optim(c(1, 1), .objf, y = .df$y)"
            )
          } else if (parms$distType == "6") {
            command <- paste0(
              "# maximum likelihood estimator\n",
              ".objf <- function(p, y) {",
              "-sum(dweibull(y, shape = p[1], scale = p[2], log = TRUE))",
              "}\n",
              ".est <- optim(c(1, 1), .objf, y = .df$y)"
            )
          }

          if (any(parms$distType == c("5", "6"))) {
            commandDoIt(command)
            registRmlist(.objf)
            registRmlist(.est)
          } else if (any(parms$distType == c("1", "2", "3", "4"))) {
            commandDoIt(command)
            registRmlist(.est)
          }

          if (any(parms$distType == c("5", "6")) && (is.na(match(".est", ls(envir = .GlobalEnv, all.names = TRUE))) || .est$convergence != 0)) {
            response <- tclvalue(RcmdrTkmessageBox(
              message = gettextKmg2("Parameter estimation failed."),
              title   = gettextKmg2("Error"),
              icon    = "error",
              type    = "ok",
              default = "ok"
            ))
            if (response == "ok") {
              errorCode <- TRUE
            }
          } else {
            .plot <- getPlot(parms)
            commandDoIt("print(.plot)")

            if (mode == 1 && parms$save == "1") savePlot(.plot)

            errorCode <- 2
          }

        }
      }
      errorCode

    },

    getGgplot = function(parms) {

      "ggplot(data = .df, aes(sample = y)) + "

    },

    getGeom = function(parms) {

      if (parms$distType == "1") {
        distParms <- "distribution = qnorm, dparams = list(mean = .est[1], sd = .est[2])"
      } else if (parms$distType == "2") {
        distParms <- "distribution = qlnorm, dparams = list(meanlog = .est[1], sdlog = .est[2])"
      } else if (parms$distType == "3") {
        distParms <- "distribution = qbeta, dparams = list(shape1 = .est[1], shape2 = .est[2])"
      } else if (parms$distType == "4") {
        distParms <- "distribution = qexp, dparams = list(rate = .est)"
      } else if (parms$distType == "5") {
        distParms <- "distribution = qgamma, dparams = list(shape = .est$par[1], scale = .est$par[2])"
      } else if (parms$distType == "6") {
        distParms <- "distribution = qweibull, dparams = list(shape = .est$par[1], scale = .est$par[2])"
      } else {
        distParms <- parms$distParms
      }

      paste0("stat_qq(", distParms, ") + ")

    },

    getMain = function(parms) {

      if (nchar(parms$main) == 0) {
        main <- ""
      } else if (parms$main == "<auto>") {
        if (parms$distType == "1") {
          main <- paste0(
            "Theoretical: qnorm(mean = ", round(.est[1], 1),
            ", sd = ", round(.est[2], 1), ")"
          )
        } else if (parms$distType == "2") {
          main <- paste0(
            "Theoretical: qlnorm(meanlog = ", round(.est[1], 1),
            ", sdlog = ", round(.est[2], 1), ")"
          )
        } else if (parms$distType == "3") {
          main <- paste0(
            "Theoretical: qbeta(shape1 = ", round(.est[1], 1),
            ", shape2 = ", round(.est[2], 1), ")"
          )
        } else if (parms$distType == "4") {
          main <- paste0("Theoretical: qexp(rate = ", round(.est, 1), ")")
        } else if (parms$distType == "5") {
          main <- paste0(
            "Theoretical: qgamma(shape = ", round(.est$par[1], 1),
            ", scale = ", round(.est$par[2], 1), ")"
          )
        } else if (parms$distType == "6") {
          main <- paste0(
            "Theoretical: qweibull(shape = ", round(.est$par[1], 1),
            ", scale = ", round(.est$par[2], 1), ")"
          )
        } else {
          main <- paste0("Theoretical: ", parms$distParms)
        }
        main <- paste0("labs(title = \"", main, "\") + ")
      } else {
        main <- paste0("labs(title = \"", parms$main, "\") + ")
      }
      main

    }

  )
)



#' Wrapper Function of Q-Q Plot Subclass
#'
#' \code{windowQQ} function is a wrapper function of \code{gqq} class for the R-commander menu bar.
#'
#' @rdname plot-gqq-windowQQ
#' @keywords hplot
#' @export
windowQQ <- function() {

  QQ <- RcmdrPlugin.KMggplot2::gqq$new()
  QQ$plotWindow()

}
