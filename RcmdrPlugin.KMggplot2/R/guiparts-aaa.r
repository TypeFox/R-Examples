#' Base Class for GUI Frames
#'
#' \code{gparts_base} class is the base class for GUI frames.
#'
#' This class is a framework for implementing subclasses which make GUI frames.
#'
#' @section Fields:
#' \describe{
#' \item{\code{frame}: }{\code{tkwin} class object; parent of widget window.} 
#' \item{\code{length}: }{Integer; number of grids.} 
#' \item{\code{back_list}: }{List of \code{tkwin} class object; list of grids.} 
#' }
#' @section Contains:
#' NULL
#' @section Methods:
#' \describe{
#' \item{\code{back(perline = 3)}: }{\code{back} method for \code{gparts_base} class.}
#' }
#' @family guiparts
#'
#' @name gparts_base-class
#' @aliases gparts_base
#' @rdname guiparts-gparts_base
#' @docType class
#' @keywords hplot
#' @export gparts_base
gparts_base <- setRefClass(

  Class = "gparts_base",
  
  fields = c("frame", "length", "back_list"),

  methods = list(

    back = function(perline = 4) {

      f <- 1
      if (!is.null(back_list)) {
        for (i in 1:length(back_list)) {
          if (i == length(back_list)) {
            do.call(
              tkgrid,
              c(back_list[f:i], list(sticky = "nw"))
            )
            do.call(
              tkgrid.configure,
              c(back_list[f:i], list(padx = c(0, 5)))
            )
          } else if (i %% (perline * 2) == 0) {
            do.call(
              tkgrid,
              c(back_list[(i - (perline * 2) + 1):i], list(sticky = "nw"))
            )
            do.call(
              tkgrid.configure,
              c(back_list[(i - (perline * 2) + 1):i], list(padx = c(0, 5)))
            )
            f <- i + 1
          }
        }
      }

      tkgrid(frame, sticky = "nw")
      tkgrid.configure(frame, pady = c(0, 5))

      return()

    }

  )
)

#' The \code{back} Method for \code{gparts_base} Class
#'
#' \code{back} method places GUI grids for R-Commander windows.
#'
#' @usage \S4method{back}{gparts_base}(perline = 3)
#' @param perline Integer; controls the number of per-line to GUI grids.
#' @family guiparts
#'
#' @name back,gparts_base-method
#' @rdname guiparts-gparts_base-back
#' @docType methods
#' @keywords hplot
NULL
