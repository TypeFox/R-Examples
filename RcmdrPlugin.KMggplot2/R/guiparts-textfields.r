#' Multiple Text Fields Subclass
#'
#' \code{textfields} class is a subclass for GUI multiple text fileds frame.
#'
#' This class is a subclass which make GUI multiple text fileds frame.
#'
#' @section Fields:
#' \describe{
#' \item{\code{frame}: }{\code{tkwin} class object; parent of widget window.}
#' \item{\code{length}: }{Integer; number of grids.}
#' \item{\code{back_list}: }{List of \code{tkwin} class object; list of grids.}
#' \item{\code{fields}: }{List of \code{textfield} class object; list of \code{textfield}.}
#' }
#' @section Contains:
#' \code{gparts_base}
#' @section Methods:
#' \describe{
#' \item{\code{back(perline = 3)}: }{\code{back} method for \code{gparts_base} class.}
#' \item{\code{front(top, initValues, titles, boxwidth = "20")}: }{\code{front} method for \code{textfields} subclass.}
#' }
#' @family guiparts
#'
#' @name textfields-class
#' @aliases textfields
#' @rdname guiparts-textfields
#' @docType class
#' @keywords hplot
#' @export textfields
textfields <- setRefClass(

  Class = "textfields",

  fields = c("fields"),

  contains = c("gparts_base"),

  methods = list(

    front = function(top, initValues, titles, boxwidth = "20") {

      if (length(initValues) != length(titles)) {
        error("length(initValues) != length(titles)")
      } else {
        length <<- length(initValues)
      }
      frame <<- tkframe(top)

      fields <<- as.list(NULL)
      for (i in 1:length) {
        fields[[i]] <<- textfield$new()
        fields[[i]]$front(frame, initValues[[i]], boxwidth, titles[[i]])
      }

      back_list <<- NULL
      for (i in 1:length) {
        back_list <<- c(back_list, list(fields[[i]]$frame))
      }

      return()

    }

  )
)



#' The \code{front} Method for \code{textfields} Subclass
#'
#' @usage \S4method{front}{textfields}(top, initValues, titles, boxwidth = "20")
#' @param top \code{tkwin} class object; top of widget window.
#' @param initValues List of character; initialization values of text fields.
#' @param titles List of character; titles of text fields.
#' @param boxwidth Character; size of text fields.
#' @family guiparts
#'
#' @name front,textfields-method
#' @rdname guiparts-textfields-front
#' @docType methods
#' @keywords hplot
NULL
