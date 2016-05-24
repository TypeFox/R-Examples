#' Radio Boxes Subclass
#'
#' \code{radioboxes} class is a subclass for GUI radio box frame.
#'
#' This class is a subclass which make GUI radio box frame.
#'
#' @section Fields:
#' \describe{
#' \item{\code{frame}: }{\code{tkwin} class object; parent of widget window.} 
#' \item{\code{length}: }{Integer; number of grids.} 
#' \item{\code{back_list}: }{List of \code{tkwin} class object; list of grids.} 
#' \item{\code{value}: }{\code{tclVar} class object; the value of the radio box frame.} 
#' \item{\code{rbradioboxes}: }{\code{tkwin} class object; the radio box frame.} 
#' }
#' @section Contains:
#' \code{gparts_base}
#' @section Methods:
#' \describe{
#' \item{\code{back(perline = 3)}: }{\code{back} method for \code{gparts_base} class.}
#' \item{\code{front(top, labels, title = "", initValue = 1, right.buttons = FALSE)}: }{\code{front} method for \code{radioboxes} subclass.}
#' }
#' @family guiparts
#'
#' @name radioboxes-class
#' @aliases radioboxes
#' @rdname guiparts-radioboxes
#' @docType class
#' @keywords hplot
#' @export radioboxes
radioboxes <- setRefClass(

  Class = "radioboxes",

  fields = c("value", "rbradioboxes"),

  contains = c("gparts_base"),

  methods = list(

    front = function(top, labels, title = "", initValue = 1, right.buttons = FALSE) {

      if (length(labels) < initValue) {
        error("length(labels) < initValue")
      } else {
        length <<- length(labels)
      }
      initValues <- 1:length
      frame <<- tkframe(top)
      value <<- tclVar(initValue)

      if (title != "") {
        tkgrid(labelRcmdr(frame, text = title, foreground = "blue"), columnspan = 2, sticky = "w")
      }
      rbradioboxes <<- as.list(NULL)
      for (i in 1:length) {
        rbradioboxes[[i]] <<- ttkradiobutton(frame, variable = value, value = initValues[i])
        if (right.buttons) {
          tkgrid(labelRcmdr(frame, text = labels[[i]], justify = "left"), rbradioboxes[[i]], sticky = "w")
        } else {
          tkgrid(rbradioboxes[[i]], labelRcmdr(frame, text = labels[[i]], justify = "left"), sticky = "w")
        }
      }

      back_list <<- list(frame)

      return()

    }

  )
)



#' The \code{front} Method for \code{radioboxes} Subclass
#'
#' @usage \S4method{front}{radioboxes}(top, labels, title = "", initValue = 1, right.buttons = FALSE)
#' @param top \code{tkwin} class object; top of widget window.
#' @param labels List of character; values of radio boxes labels.
#' @param title Character; the title of the radio box frame.
#' @param initValue Integer; the initialization value of the radio box frame.
#' @param right.buttons Boolean; whether radio boxes are right aligned or not.
#' @family guiparts
#'
#' @name front,radioboxes-method
#' @rdname guiparts-radioboxes-front
#' @docType methods
#' @keywords hplot
NULL
