#' Check Boxes Subclass
#'
#' \code{checkboxes} class is a subclass for GUI check box frame.
#'
#' This class is a subclass which make GUI check box frame.
#'
#' @section Fields:
#' \describe{
#' \item{\code{frame}: }{\code{tkwin} class object; parent of widget window.} 
#' \item{\code{length}: }{Integer; number of grids.} 
#' \item{\code{back_list}: }{List of \code{tkwin} class object; list of grids.} 
#' \item{\code{value}: }{List of \code{tclVar} class object; values of the check box frame.} 
#' \item{\code{cbcheckboxes}: }{\code{tkwin} class object; the check box frame.} 
#' }
#' @section Contains:
#' \code{gparts_base}
#' @section Methods:
#' \describe{
#' \item{\code{back(perline = 3)}: }{\code{back} method for \code{gparts_base} class.}
#' \item{\code{front(top, initValues, labels, title = "", right.buttons = FALSE)}: }{\code{front} method for \code{checkboxes} subclass.}
#' }
#' @family guiparts
#'
#' @name checkboxes-class
#' @aliases checkboxes
#' @rdname guiparts-checkboxes
#' @docType class
#' @keywords hplot
#' @export checkboxes
checkboxes <- setRefClass(

  Class = "checkboxes",

  fields = c("value", "cbcheckboxes"),

  contains = c("gparts_base"),

  methods = list(

    front = function(top, initValues, labels, title = "", right.buttons = FALSE) {

      if (length(initValues) != length(labels)) {
        error("length(initValues) != length(labels)")
      } else {
        length <<- length(initValues)
      }
      frame <<- tkframe(top)

      if (title != "") {
        tkgrid(labelRcmdr(frame, text = title, foreground = "blue"), columnspan = 2, sticky = "w")
      }
      value <<- as.list(NULL)
      cbcheckboxes <<- as.list(NULL)
      for (i in 1:length) {
        value[[i]] <<- tclVar(initValues[[i]])
        cbcheckboxes[[i]] <<- tkcheckbutton(frame, variable = value[[i]])
        if (right.buttons) {
          tkgrid(labelRcmdr(frame, text = labels[[i]], justify = "left"), cbcheckboxes[[i]], sticky = "w")
        } else {
          tkgrid(cbcheckboxes[[i]], labelRcmdr(frame, text = labels[[i]], justify = "left"), sticky = "w")
        }
      }

      back_list <<- list(frame)

      return()

    }

  )
)



#' The \code{front} Method for \code{checkboxes} Subclass
#'
#' @usage \S4method{front}{checkboxes}(top, initValues, labels, title = "", right.buttons = FALSE)
#' @param top \code{tkwin} class object; top of widget window.
#' @param initValues List of integer; initialization values of check boxes.
#' @param labels List of character; values of check boxes labels.
#' @param title Character; the title of the check box frame.
#' @param right.buttons Boolean; whether check boxes are right aligned or not.
#' @family guiparts
#'
#' @name front,checkboxes-method
#' @rdname guiparts-checkboxes-front
#' @docType methods
#' @keywords hplot
NULL
