#' Factorize Subclass
#'
#' \code{factorize} class is a subclass for factorizing numeric variables.
#'
#' This class is a subclass which show dialog boxes of a factorizer for graphics editing.
#'
#' @section Fields:
#' \describe{
#' \item{\code{top}: }{\code{tkwin} class object; parent of widget window.}
#' \item{\code{alternateFrame}: }{\code{tkwin} class object; a special frame for some GUI parts.}
#' \item{\code{vbbox1}: }{\code{variableboxes} class object; the frame to select variables.}
#' \item{\code{rbbox1}: }{\code{radioboxes} class object; the frame to set options.}
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
#' @export factorize
#' @name factorize-class
#' @aliases factorize
#' @rdname plot-factorize
#' @docType class
#' @keywords hplot
factorize <- setRefClass(
  
  Class = "factorize",
  
  fields = c("vbbox1", "rbbox1"),
  
  contains = c("plot_base"),
  
  methods = list(
    
    #' Plot Windows
    plotWindow = function() {
      
      # note: The initializeDialog() generates "top"
      initializeDialog(window = topwindow, title = getWindowTitle())
      top            <<- topwindow
      alternateFrame <<- tkframe(top)
      
      setFront()
      
      parms <- getParms()
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

        }
        
        removeRmlist()
        
        # tkinsert(LogWindow(), "end", codes)
        
        activateMenus()
        tkfocus(CommanderWindow())
        
      }
      
      setBack()
      
      # note: The OKCancelHelp() generates "buttonsFrame"
      OKCancelHelp(window = top, helpSubject = getHelp())
      
      tkgrid(buttonsFrame, sticky = "nw")
      dialogSuffix()
      
      return()
      
    },
    
    setFront = function() {
      
      vbbox1 <<- variableboxes$new()
      vbbox1$front(
        top       = top, 
        types     = list(Numeric()),
        titles    = list(
          gettextKmg2("Variable (pick one or more)")
        ),
        modes     = list("multiple"),
        initialSelection = list(0)
      )
      
      rbbox1 <<- radioboxes$new()
      rbbox1$front(
        top    = top,
        labels = list(
          gettextKmg2("Simple factorization"),
          gettextKmg2("Categorize according to the quartiles")
        ),
        title  = gettextKmg2("Options")
      )
      
    },
    
    setBack = function() {
      
      vbbox1$back()
      rbbox1$back()
      
    },
    
    getWindowTitle = function() {
      
      gettextKmg2("Factorizing numeric variables")
      
    },
    
    getHelp = function() {
      
      "factor"
      
    },
    
    getParms = function() {
      
      x      <- getSelection(vbbox1$variable[[1]])
      y      <- character(0)
      z      <- character(0)
      
      s      <- character(0)
      t      <- character(0)
      
      xlab   <- character(0)
      xauto  <- character(0)
      ylab   <- character(0)
      yauto  <- character(0)
      zlab   <- character(0)
      main   <- character(0)
      
      size   <- character(0)
      family <- character(0)
      colour <- character(0)
      save   <- character(0)
      theme  <- character(0)
      
      factoriseType <- tclvalue(rbbox1$value)
      
      list(
        x = x, y = y, z = z, s = s, t = t,
        xlab = xlab, xauto = xauto, ylab = ylab, yauto = yauto, zlab = zlab, main = main,
        size = size, family = family, colour = colour, save = save, theme = theme,
        factoriseType = factoriseType
      )
      
    },
    
    checkError = function(parms) {
      
      if (length(parms$x) == 0) {
        errorCondition(
          recall  = windowScatter,
          message = gettextKmg2("Variables are not selected")
        )
        errorCode <- TRUE
      } else {
        setDataframe(parms)
        
        errorCode <- 2
      }
      errorCode
      
    },
    
    setDataframe = function(parms) {
      
      if (parms$factoriseType == "1") {
        command <- paste0(
          ActiveDataSet(), "$", parms$x, ".f",
          " <- ",
          "factor(", ActiveDataSet(), "$", parms$x, ")"
        )
      } else {
        command <- paste0(
          ActiveDataSet(), "$", parms$x, ".fq",
          " <- ",
          "factor(cut(", ActiveDataSet(), "$", parms$x, ", breaks = quantile(", ActiveDataSet(), "$", parms$x, ", c(0, 0.25, 0.50, 0.75, 1)), include.lowest = TRUE))"
        )
      }
      
      doItAndPrint(command)
      activeDataSet(ActiveDataSet())
      
    }
    
  )
)



#' Wrapper Function of Factorize Subclass
#'
#' \code{windowFactorize} function is a wrapper function of \code{factorize} class for the R-commander menu bar.
#'
#' @rdname plot-factorize-windowFactorize
#' @keywords hplot
#' @export
windowFactorize <- function() {
  
  Factorize <- RcmdrPlugin.KMggplot2::factorize$new()
  Factorize$plotWindow()
  
}
