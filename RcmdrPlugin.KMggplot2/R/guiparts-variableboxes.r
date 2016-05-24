#' Variable Boxes Subclass
#'
#' \code{variableboxes} class is a subclass for GUI variable box frame.
#'
#' This class is a subclass which make GUI variable box frame.
#'
#' @section Fields:
#' \describe{
#' \item{\code{frame}: }{\code{tkwin} class object; parent of widget window.} 
#' \item{\code{length}: }{Integer; number of grids.} 
#' \item{\code{back_list}: }{List of \code{tkwin} class object; list of grids.} 
#' \item{\code{variable}: }{List of \code{listbox} class object; variables of the variable box frame.} 
#' }
#' @section Contains:
#' \code{gparts_base}
#' @section Methods:
#' \describe{
#' \item{\code{back(perline = 3)}: }{\code{back} method for \code{gparts_base} class.}
#' \item{\code{front(top, types, titles, modes = "default")}: }{\code{front} method for \code{variableboxes} subclass.}
#' }
#' @family guiparts
#'
#' @name variableboxes-class
#' @aliases variableboxes
#' @rdname guiparts-variableboxes
#' @docType class
#' @keywords hplot
#' @export variableboxes
variableboxes <- setRefClass(

  Class = "variableboxes",

  fields = c("variable"),

  contains = c("gparts_base"),

  methods = list(

    front = function(top, types, titles, initialSelection = FALSE, modes = "default") {

      if (length(types) != length(titles)) {
        error("length(types) != length(titles)")
      } else {
        length <<- length(types)
      }
      if (length(modes) == 1) {
        if (modes == "default") {
          modes <- lapply(1:length, function(x) "single")
        }    
      }
      if (length(types) != length(modes)) {
        error("length(types) != length(modes)")
      }   
      if (length(initialSelection) == 1) {
        initialSelection <- lapply(1:length, function(x) initialSelection)
      } else if (length(initialSelection) != length) {
        error("length(initialSelection) != length(types)")
      }
      frame <<- tkframe(top)

      variable <<- lapply(1:length, function(i, types, modes, initialSelection, titles) {
          variableListBox(
            frame,
            variableList     = unlist(types[i]),
            selectmode       = unlist(modes[i]),
            initialSelection = unlist(initialSelection[i]),
            title            = unlist(titles[i])
          )
        }, types, modes, initialSelection, titles
      )

      back_list <<- NULL
      for (i in 1:length) {
        back_list <<- c(back_list, list(variable[[i]]$frame))
      }

      return()

    }

  )
)



#' The \code{front} Method for \code{variableboxes} Subclass
#'
#' @usage \S4method{front}{variableboxes}(top, types, titles, modes = "default")
#' @param top \code{tkwin} class object; top of widget window.
#' @param types List of character; types of variableLists.
#' @param titles List of character; titles of variable box frame.
#' @param modes List of character; select modes.
#' @family guiparts
#'
#' @name front,variableboxes-method
#' @rdname guiparts-variableboxes-front
#' @docType methods
#' @keywords hplot
NULL
