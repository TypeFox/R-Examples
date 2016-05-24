#' Text Field Subclass
#'
#' \code{textfield} class is a subclass for GUI text field frame.
#'
#' This class is a subclass which make GUI text field frame.
#'
#' @section Fields:
#' \describe{
#' \item{\code{frame}: }{\code{tkwin} class object; parent of widget window.}
#' \item{\code{length}: }{Integer; number of grids.}
#' \item{\code{back_list}: }{List of \code{tkwin} class object; list of grids.}
#' \item{\code{value}: }{\code{tclVar} class object; the value of the text field frame.}
#' }
#' @section Contains:
#' \code{gparts_base}
#' @section Methods:
#' \describe{
#' \item{\code{back(perline = 3)}: }{\code{back} method for \code{gparts_base} class.}
#' \item{\code{front(top, initialValue, boxwidth, title)}: }{\code{front} method for \code{textfield} subclass.}
#' }
#' @family guiparts
#'
#' @name textfield-class
#' @aliases textfield
#' @rdname guiparts-textfield
#' @docType class
#' @keywords hplot
#' @export textfield
textfield <- setRefClass(

  Class = "textfield",

  fields = c("value"),

  contains = c("gparts_base"),

  methods = list(

    front = function(top, initialValue, boxwidth, title) {

      frame <<- tkframe(top)
      value <<- tclVar(initialValue)
      textboxEntry <- ttkentry(frame, width=boxwidth, textvariable=value)
      tkgrid(labelRcmdr(frame, text=title, fg="blue"), sticky="nw")
      tkgrid(textboxEntry, sticky="nw")

      back_list <<- as.list(NULL)

      return()

    }

  )
)



#' The \code{front} Method for \code{textfield} Subclass
#'
#' @usage \S4method{front}{textfield}(top, initialValue, boxwidth, title)
#' @param top \code{tkwin} class object; top of widget window.
#' @param initialValue Character; initialization values of the text field frame.
#' @param boxwidth Integer; the size of the text field frame.
#' @param title Character; the title of the text field frame.
#' @family guiparts
#'
#' @name front,textfield-method
#' @rdname guiparts-textfield-front
#' @docType methods
#' @keywords hplot
NULL
