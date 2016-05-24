#' Kaplan-Meier Plot Subclass
#'
#' \code{gkm} class is a subclass for Kaplan-Meier plots.
#'
#' This class is a subclass which show dialog boxes of Kaplan-Meier plots for graphics editing.
#'
#' @section Fields:
#' \describe{
#' \item{\code{top}: }{\code{tkwin} class object; parent of widget window.}
#' \item{\code{alternateFrame}: }{\code{tkwin} class object; a special frame for some GUI parts.}
#' \item{\code{vbbox1}: }{\code{variableboxes} class object; the frame to select variables.}
#' \item{\code{vbbox2}: }{\code{variableboxes} class object; the frame to select facet variables.}
#' \item{\code{lbbox1}: }{\code{textfields} class object; the frame to set axis labels, the legend label, and the main title.}
#' \item{\code{lbbox2}: }{\code{textfields} class object; the frame to set the tick count.}
#' \item{\code{rbbox1}: }{\code{radioboxes} class object; the frame to set the plot type.}
#' \item{\code{cbbox1}: }{\code{checkboxes} class object; the frame to set option(s).}
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
#' @name gkm-class
#' @aliases gkm
#' @rdname plot-gkm
#' @docType class
#' @keywords hplot
#' @importFrom grid grid.newpage pushViewport viewport
#' @importFrom survival survfit Surv survdiff
#' @importFrom plyr ddply .
#' @export gkm
gkm <- setRefClass(

  Class = "gkm",

  fields = c("vbbox1", "vbbox2", "lbbox1", "lbbox2", "rbbox1", "cbbox1", "tbbox1"),

  contains = c("plot_base"),

  methods = list(

    setFront = function() {

      vbbox1 <<- variableboxes$new()
      vbbox1$front(
        top       = top, 
        types     = list(nonFactors(), Variables(), Factors()),
        titles    = list(
          gettextKmg2("Time variable"),
          gettextKmg2("Event variable (0=censor, 1=event)"),
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
        initValues = list("Time from entry", "Proportion of survival", "<auto>", ""),
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
          gettextKmg2("None"),
          gettextKmg2("Inside"),
          gettextKmg2("Outside (facet ignored)")
        ),
        title  = gettextKmg2("No at risk.")
      )

      lbbox2 <<- textfields$new()
      lbbox2$front(
        top        = alternateFrame,
        initValues = list("4"),
        titles     = list(
          gettextKmg2("Tick count")
        )
      )

      cbbox1 <<- checkboxes$new()
      cbbox1$front(
        top        = alternateFrame,
        initValues = list("0", "0", "0", "0", "0"),
        labels     = list(
          gettextKmg2("Confidence interval (dotted lines)"),
          gettextKmg2("Confidence interval (band)"),
          gettextKmg2("Dot censored symbol"),
          gettextKmg2("P-value (log-rank test)"),
          gettextKmg2("Reference line at median survival")
        ),
        title      = gettextKmg2("Options")
      )

      tbbox1 <<- toolbox$new()
      tbbox1$front(top)

    },

    setBack = function() {

      vbbox1$back()
      vbbox2$back()
      lbbox1$back(4)

      boxlist <- c(
        list(rbbox1$frame),
        list(labelRcmdr(alternateFrame, text="    ")),
        list(lbbox2$frame),
        list(labelRcmdr(alternateFrame, text="    ")),
        list(cbbox1$frame)
      )
      do.call(tkgrid, c(lbbox2$back_list, list(sticky="nw")))
      do.call(tkgrid, c(boxlist, list(sticky="nw")))
      tkgrid(alternateFrame, stick="nw")
      tkgrid(labelRcmdr(alternateFrame, text="    "), stick="nw")

      tbbox1$back()

    },

    getWindowTitle = function() {
      
      gettextKmg2("Kaplan-Meier plot")
      
    },
    
    getHelp = function() {
      
      "survfit"
      
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
      xauto  <- "Time from entry"
      ylab   <- tclvalue(lbbox1$fields[[2]]$value)
      yauto  <- "Proportion of survival"
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
      
      plotType  <- tclvalue(rbbox1$value)
      tickCount <- tclvalue(lbbox2$fields[[1]]$value)
      confInt   <- tclvalue(cbbox1$value[[1]])
      confIntB  <- tclvalue(cbbox1$value[[2]])
      dotCensor <- tclvalue(cbbox1$value[[3]])
      pValue    <- tclvalue(cbbox1$value[[4]])
      refMedian <- tclvalue(cbbox1$value[[5]])
      
      if (is.na(as.numeric(tickCount)) || as.numeric(tickCount) <= 0) {
        tickCount <- "4"
      }
      if (plotType == "3") {
        s <- character(0)
        t <- character(0)
      }

      if (length(s) != 0 && length(t) != 0) {
        zst <- "\"z\", \"s\", \"t\""
      } else if (length(s) != 0) {
        zst <- "\"z\", \"s\""
      } else if (length(t) != 0) {
        zst <- "\"z\", \"t\""
      } else {
        zst <- "\"z\""
      }

      list(
        x = x, y = y, z = z, s = s, t = t,
        xlab = xlab, xauto = xauto, ylab = ylab, yauto = yauto, zlab = zlab, main = main,
        size = size, family = family, colour = colour, save = save, theme = theme,
        plotType = plotType, tickCount = tickCount, confInt = confInt, confIntB = confIntB,
        dotCensor = dotCensor, pValue = pValue, refMedian = refMedian, zst = zst
      )

    },

    checkError = function(parms) {

      if (length(parms$x) == 0) {
        errorCondition(
          recall  = windowScatter,
          message = gettextKmg2("Time variable is not selected.")
        )
        errorCode <- TRUE
      } else if (length(parms$y) == 0) {
        errorCondition(
          recall  = windowScatter,
          message = gettextKmg2("Event variable is not selected.")
        )
        errorCode <- TRUE
      } else {
        
        if (mode == 1) {
          logger("require(\"ggplot2\")")
        }

        setDataframe(parms)

        strataNum <- length(unique(.df$z))

        if (length(parms$s) != 0 && length(parms$t) != 0) {
          command <- ".fit <- survival::survfit(survival::Surv(time = x, event = y, type = \"right\") ~ z + s + t, .df)"
        } else if (length(parms$s) != 0) {
          command <- ".fit <- survival::survfit(survival::Surv(time = x, event = y, type = \"right\") ~ z + s, .df)"
        } else if (length(parms$t) != 0) {
          command <- ".fit <- survival::survfit(survival::Surv(time = x, event = y, type = \"right\") ~ z + t, .df)"
        } else {
          command <- ".fit <- survival::survfit(survival::Surv(time = x, event = y, type = \"right\") ~ z, .df)"
        }
        response <- tryCatch({
            commandDoIt(command)
            ""
          }, error = function(ex) {
            tclvalue(RcmdrTkmessageBox(
              message = gettextKmg2("Estimation failed in survfit().  This data may be inappropriate."),
              title   = gettextKmg2("Error"),
              icon    = "error",
              type    = "ok",
              default = "ok"
            ))
          }
        )
        if (response == "ok") {
          return(TRUE)
        }

        registRmlist(.fit)
        
        
        if (parms$pValue == "1" && strataNum > 1) {
          
          if (length(parms$s) != 0 && length(parms$t) != 0) {
            command <- "s, t"
          } else if (length(parms$s) != 0) {
            command <- "s"
          } else if (length(parms$t) != 0) {
            command <- "t"
          } else {
            command <- ""
          }
          
          if (parms$plotType == "2") {
            ypos <- 0.25
          } else {
            ypos <- 0
          }
          
          command <- paste0(
            ".pval <- plyr::ddply(.df, plyr::.(", command, "),\n",
            " function(x) {\n",
            "  data.frame(\n",
            "   x = 0, y = ", ypos,", df = ", strataNum - 1, ",\n",
            "   chisq = survival::survdiff(\n",
            "    survival::Surv(time = x, event = y, type = \"right\") ~ z, x\n",
            "   )$chisq\n",
            ")})\n",
            ".pval$label <- paste0(\n",
            "  \"paste(italic(p), \\\" = \",\n",
            "  signif(1 - pchisq(.pval$chisq, .pval$df), 3),\n",
            "  \"\\\")\"\n",
            ")"
            )
          registRmlist(.pval)
          commandDoIt(command)

          
          if (length(parms$s) != 0 && length(parms$t) != 0) {
            command <- "\"s\", \"t\""
          } else if (length(parms$s) != 0) {
            command <- "\"s\""
          } else if (length(parms$t) != 0) {
            command <- "\"t\""
          } else {
            command <- ""
          }
          if (command != "") {
            command <- paste0(
              ".ppval <- by(.df, .df[, c(", command, ")], \n",
              " function(x) {\n",
              "   survival::survdiff(\n",
              "    survival::Surv(time = x, event = y, type = \"right\") ~ z, x\n",
              "   )\n",
              "})\n"
            )
            registRmlist(.ppval)
          } else {
            command <- ".ppval <- survival::survdiff(survival::Surv(time = x, event = y, type = \"right\") ~ z, .df)"
            registRmlist(.ppval)
          }
          justDoIt(command)
          
        } else if (parms$pValue == "1" && strataNum == 1) {
          
          Message(message = gettextKmg2("There is only 1 group. P-value could not be calculated."),
                  type = "warning")
          
        }
        
        if (parms$refMedian == "1") {
          
          if (length(parms$s) != 0 && length(parms$t) != 0) {
            command <- ".pmed <- survival::survfit(survival::Surv(time = x, event = y, type = \"right\") ~ z + s + t, .df)"
          } else if (length(parms$s) != 0) {
            command <- ".pmed <- survival::survfit(survival::Surv(time = x, event = y, type = \"right\") ~ z + s, .df)"
          } else if (length(parms$t) != 0) {
            command <- ".pmed <- survival::survfit(survival::Surv(time = x, event = y, type = \"right\") ~ z + t, .df)"
          } else {
            command <- ".pmed <- survival::survfit(survival::Surv(time = x, event = y, type = \"right\") ~ z, .df)"
          }
          registRmlist(.pmed)
          justDoIt(command)
          
        }
        
        if (length(parms$s) != 0 && length(parms$t) != 0) {
          command <- "\"x\", \"z\", \"s\", \"t\""
        } else if (length(parms$s) != 0) {
          command <- "\"x\", \"z\", \"s\""
        } else if (length(parms$t) != 0) {
          command <- "\"x\", \"z\", \"t\""
        } else {
          command <- "\"x\", \"z\""
        }
        command <- paste0(
          ".fit <- data.frame(",
            "x = .fit$time, y = .fit$surv, nrisk = .fit$n.risk, ",
            "nevent = .fit$n.event, ncensor= .fit$n.censor, ",
            "upper = .fit$upper, lower = .fit$lower",
          ")\n",
          ".df <- .df[!duplicated(.df[,c(", command, ")]), ]\n",
          ".df <- .fit <- data.frame(.fit, .df[, c(", parms$zst, "), drop = FALSE])"
        )
        commandDoIt(command)
        
        if (parms$refMedian == "1") {
          
          if (length(parms$s) != 0 && length(parms$t) != 0) {
            command <- "z, s, t"
          } else if (length(parms$s) != 0) {
            command <- "z, s"
          } else if (length(parms$t) != 0) {
            command <- "z, t"
          } else {
            command <- "z"
          }
          
          command <- paste0(
            ".med <- plyr::ddply(.fit, plyr::.(", command, "), function(x) {\n",
            "data.frame(\n",
            " median = min(subset(x, y < (0.5 + .Machine$double.eps^0.5))$x)\n",
            ")})"
          )
          
          registRmlist(.med)
          commandDoIt(command)
          
        }

        if (!all(round(.fit$x - .fit$x.1, 1e-13) == 0)) {
          tclvalue(RcmdrTkmessageBox(
            message = gettextKmg2("Unknown error.  This plot may be broken."),
            title   = gettextKmg2("Error"),
            icon    = "error",
            type    = "ok",
            default = "ok"
          ))
        }

        .plot <- getPlot(parms)
        if (mode == 1) {
          logger("print(.plot)")
        }
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
          return(TRUE)
        }

        if (mode == 1 && parms$save == "1") savePlot(.plot)

        errorCode <- 2
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
      } else {
        var <- c(var, paste0("z = factor(\"At risk\")"))
      }
      if (length(parms$s) != 0) {
        var <- c(var, paste0("s = ", ActiveDataSet(), "$", parms$s))
      }
      if (length(parms$t) != 0) {
        var <- c(var, paste0("t = ", ActiveDataSet(), "$", parms$t))
      }
      command <- do.call(paste, c(var, list(sep = ", ")))
      command <- paste0(".df <- na.omit(data.frame(", command, "))")
      commandDoIt(command)

      command <- paste0(
        ".df <- .df[do.call(order, .df[, c(", parms$zst, ", \"x\"), drop = FALSE]), , drop = FALSE]"
      )
      commandDoIt(command)

      registRmlist(.df)

    },

    getGgplot = function(parms) {

      "ggplot(data = .fit, aes(x = x, y = y, colour = z)) + "

    },

    getGeom = function(parms) {

      if (parms$confInt == "1") {
        geom <- paste0(
          "geom_step(data = subset(.fit, !is.na(upper)), aes(y = upper), size = 1, lty = 2, alpha = 0.5, show.legend = FALSE, na.rm = FALSE) + ",
          "geom_step(data = subset(.fit, !is.na(lower)), aes(y = lower), size = 1, lty = 2, alpha = 0.5, show.legend = FALSE, na.rm = FALSE) + "
        )
      } else {
        geom <- ""
      }
      
      if (parms$confIntB == "1") {
        geom <- paste0(
          geom,
          "RcmdrPlugin.KMggplot2::geom_stepribbon(data = .fit, aes(x = x, ymin = lower, ymax = upper, fill = z), alpha = 0.25, colour = \"transparent\", show.legend = FALSE, kmplot = TRUE) + "
        )
      }
      
      geom <- paste0(geom, "geom_step(size = 1.5) + ")

      if (nrow(.cens) > 0) {
        if (parms$dotCensor == "1") {
          geom <- paste0(
            geom,
            "geom_point(data = .cens, aes(x = x, y = y, colour = z), size = 3.5) + ",
            "geom_point(data = .cens, aes(x = x, y = y), size = 2, color = \"white\") + "
          )
        } else {
          geom <- paste0(
            geom,
            "geom_linerange(data = .cens, aes(x = x, ymin = y, ymax = y + 0.02), size = 1.5) + "
          )
        }
      }

      if (parms$plotType == "2") {
        geom <- paste0(
          geom,
          "geom_text(data = .nrisk, aes(y = y, x = x, label = Freq, colour = z), show.legend = FALSE, size = ", parms$size , " * 0.282, family = \"", parms$family, "\") + "
        )
      }
       
      if (parms$pValue == "1" && length(unique(.df$z)) > 1) {
        geom <- paste0(
          geom,
          "geom_text(data = .pval, aes(y = y, x = x, label = label), colour = \"black\", hjust = 0, vjust = -0.5, parse = TRUE, show.legend = FALSE, size = ", parms$size , " * 0.282, family = \"", parms$family, "\") + "
        )
      }
      
      if (parms$refMedian == "1" && !all(is.infinite(.med$median))) {
        geom <- paste0(
          geom,
          "geom_vline(data = .med, aes(xintercept = median), colour = \"black\", lty = 2) + "
        )
      } else if (parms$refMedian == "1" && all(is.infinite(.med$median))) {
        Message(message = gettextKmg2("Median survival times did not exist."),
                type = "warning")
      }
      
      geom

    },

    getScale = function(parms) {
      
      nTick    <- as.numeric(parms$tickCount) - 1
      byLength <- signif((max(.fit$x) + nTick/2)/nTick, -round(log10(max(.fit$x)), 0))
      scale <- paste0(
        "scale_x_continuous(breaks = seq(0, ", byLength * nTick, ", by = ", byLength, "), limits = c(0, ", byLength * nTick, ")) + ", 
        "scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0)) + ",
        "scale_colour_brewer(palette = \"", parms$colour, "\") + "
      )
      if (parms$confIntB == "1") {
        scale <- paste0(
          scale,
          "scale_fill_brewer(palette = \"", parms$colour, "\") + "
        )
      }
      scale
      
    },

    getZlab = function(parms) {

      if (length(parms$z) == 0) {
        zlab <- ""
      } else if (parms$zlab == "<auto>") {
        zlab <- paste0("labs(colour = \"", parms$z, "\") + ")
      } else {
        zlab <- paste0("labs(colour = \"", parms$zlab, "\") + ")
      }
      zlab

    },

    getOpts = function(parms) {

      opts <- list()
      if (length(parms$s) != 0 || length(parms$t) != 0) {
        opts <- c(opts, "panel.margin = unit(0.3, \"lines\")")
      }

      if (length(parms$z) == 0) {
        opts <- c(opts, "legend.position = \"none\"")
      } else {
        if (nchar(parms$zlab) == 0) {
          opts <- c(opts, "legend.title = element_blank()")
        }
        if (parms$plotType == "3") {
          opts <- c(opts, "legend.position = c(1, 1)", "legend.justification = c(1, 1)")
        } else {
          opts <- c(opts, "legend.position = \"right\"")
        }
      }

      if (length(opts) != 0) {
        opts <- do.call(paste, c(opts, list(sep = ", ")))
        opts <- paste0(" + theme(", opts, ")")
      } else {
        opts <- ""
      }
      opts

    },

    getPlot = function(parms) {

      command <- paste0(
        ".df <- .fit <- rbind(unique(data.frame(x = 0, y = 1, nrisk = NA, nevent = NA, ncensor = NA, upper = 1, lower = 1, .df[, c(", parms$zst, "), drop = FALSE])), .fit)"
      )
      commandDoIt(command)

      command <- ".cens <- subset(.fit, ncensor == 1)"
      commandDoIt(command)
      registRmlist(.cens)

      nTick    <- as.numeric(parms$tickCount) - 1
      byLength <- signif((max(.fit$x) + nTick/2)/nTick, -round(log10(max(.fit$x)), 0))
      
      if (parms$plotType != "1") {

        command <- paste0(
          ".tmp1 <- data.frame(as.table(by(.df, .df[, c(", parms$zst, "), drop = FALSE], function(d) max(d$nrisk, na.rm = TRUE))))\n",
          ".tmp1$x <- 0\n",
          ".nrisk <- .tmp1\n",
          "for (i in 1:", nTick, ") {",
          ".df <- subset(.fit, x < ", byLength, " * i); ",
          ".tmp2 <- data.frame(as.table(by(.df, .df[, c(", parms$zst, "), drop = FALSE], function(d) if (all(is.na(d$nrisk))) NA else min(d$nrisk - d$nevent - d$ncensor, na.rm = TRUE)))); ",
          ".tmp2$x <- ", byLength, " * i; ",
          ".tmp2$Freq[is.na(.tmp2$Freq)] <- .tmp1$Freq[is.na(.tmp2$Freq)]; ",
          ".tmp1 <- .tmp2; ",
          ".nrisk <- rbind(.nrisk, .tmp2)",
          "}"
        )
        commandDoIt(command)
        command <- paste0(
          ".nrisk$y <- rep(seq(", 0.05 * length(unique(factor(.fit$z))) - 0.025, ", 0.025, -0.05), ", nrow(.nrisk) / length(unique(factor(.fit$z))), ")"
        )
        commandDoIt(command)
        registRmlist(.tmp1)
        registRmlist(.tmp2)
        registRmlist(.nrisk)
        
      }

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

      if (parms$plotType == "3") {

        if (length(parms$z) != 0) {
          opts <- " + theme(legend.position = \"right\")"
          if (length(unique(factor(.fit$z))) == 2) {
            command <- ".nrisk$y <- ((.nrisk$y - 0.025) / (max(.nrisk$y) - 0.025) + 0.5) * 0.5"
          } else {
            command <- ".nrisk$y <- ((.nrisk$y - 0.025) / (max(.nrisk$y) - 0.025) + 0.125) * 0.8"
          }
        } else {
          command <- ".nrisk$y <- 0.5"
        }
        commandDoIt(command)
        
        command <- paste0(
          ".plot2 <- ggplot(data = .nrisk, aes(x = x, y = y, label = Freq, colour = z)) + ",
          "geom_text(size = ", parms$size , " * 0.282, family = \"", parms$family, "\") + ",
          "scale_x_continuous(breaks = seq(0, ", byLength * nTick, ", by = ", byLength, "), limits = c(0, ", byLength * nTick, ")) + ", 
          "scale_y_continuous(limits = c(0, 1)) + ",
          "scale_colour_brewer(palette = \"", parms$colour, "\") + ",
          ylab,
          "RcmdrPlugin.KMggplot2::theme_natrisk(", parms$theme, ", ", parms$size, ", \"", parms$family, "\")"
        )
        commandDoIt(command)

        command <- paste0(
          ".plot3 <- ggplot(data = subset(.nrisk, x == 0), aes(x = x, y = y, label = z, colour = z)) + ",
          "geom_text(hjust = 0, size = ", parms$size , " * 0.282, family = \"", parms$family, "\") + ",
          "scale_x_continuous(limits = c(-5, 5)) + ",
          "scale_y_continuous(limits = c(0, 1)) + ",
          "scale_colour_brewer(palette = \"", parms$colour, "\") + ",
          "RcmdrPlugin.KMggplot2::theme_natrisk21(", parms$theme, ", ", parms$size, ", \"", parms$family, "\")"
        )
        commandDoIt(command)
        
        command <- paste0(
          ".plotb <- ggplot(.df, aes(x = x, y = y)) + geom_blank() + ",
          "RcmdrPlugin.KMggplot2::theme_natriskbg(", parms$theme, ", ", parms$size, ", \"", parms$family, "\")"
        )
        commandDoIt(command)
        
        registRmlist(.plot2)
        registRmlist(.plot3)
        registRmlist(.plotb)
        
        command <- paste0(
          "grid::grid.newpage(); ",
          "grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 2, ",
          "heights = unit(c(1, ", length(unique(.fit$z))/2 + 2, "), c(\"null\", \"lines\")), ",
          "widths  = unit(c(4, 1), c(\"lines\", \"null\"))))); ",
          "print(.plotb, vp = grid::viewport(layout.pos.row = 1:2, layout.pos.col = 1:2)); ",
          "print(.plot , vp = grid::viewport(layout.pos.row = 1  , layout.pos.col = 1:2)); ",
          "print(.plot2, vp = grid::viewport(layout.pos.row = 2  , layout.pos.col = 1:2)); ",
          "print(.plot3, vp = grid::viewport(layout.pos.row = 2  , layout.pos.col = 1  )); ",
          ".plot <- recordPlot()"
        )
        commandDoIt(command)
      
      }
      
      if (parms$refMedian == "1" && !all(is.infinite(.med$median))) {
        commandDoIt("print(.pmed)", log = FALSE)
        commandDoIt("cat(\"\n\")", log = FALSE)
      }
      
      if (parms$pValue == "1" && length(unique(.df$z)) > 1) {
        commandDoIt("print(.ppval)", log = FALSE)
        commandDoIt("cat(\"\n\")", log = FALSE)
      }
      
      return(.plot)

    }

  )
)



#' Wrapper Function of Kaplan-Meier Plot Subclass
#'
#' \code{windowKM} function is a wrapper function of \code{gkm} class for the R-commander menu bar.
#'
#' @rdname plot-gkm-windowKM
#' @keywords hplot
#' @export
windowKM <- function() {

  KM <- RcmdrPlugin.KMggplot2::gkm$new()
  KM$plotWindow()

}
