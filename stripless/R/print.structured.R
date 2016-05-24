#' @title Methods for Structured Objects
#' @description Print/plot and summary methods for class 'structured' objects.
#'
#' @details
#' The \code{print} and \code{plot} methods produce a plot and informative legend
#' for "structured" objects. The \code{plot} method is an alias for the print
#' method. The \code{summary} method gives a simple summary of the object.
#'
#' @aliases plot.structured
#' @param x,object The object to be displayed or summarized.
#'
#' @param legendLoc An optional location of a legend describing the plot layout
#'   to be used as a legend on the trellis plot. Note that a text legend is
#'   \emph{always} printed on the console. If used, this argument specified the
#'   position of a legend in the trellis display. It must be one of
#'   \code{"left"}, \code{"right"}, \code{"top"}, \code{"bottom"}, or
#'   \code{"newpage"} (matching is done using \code{\link{match.arg}}, so
#'   \code{legendLoc} can be abbreviated). Any of the first four of these will
#'   become the name of a component of the \code{legend} argument of the
#'   \code{\link[lattice]{xyplot}} call, and so must not conflict with any in a
#'   \code{legend} argument that may already be part of the
#'   \code{\link{strucplot}} call.
#'
#' If \code{legendLoc = "newpage"}, the legend will be plotted centered on a new
#' trellis page.
#'
#' @param legend A \emph{function} that constructs a grob to use as a plot legend. It
#' must accept at least 2 arguments, \code{"struc"} and \code{"legendLoc"}, which
#' will be passed the \code{"structure"} attribute of the object to
#' be plotted and the \code{legendLoc} argument above.
#' Additional arguments may be passed in the \dots argument of the print
#' call (below).
#'
#' The default is the unexported \code{\link{defaultStrucLegend}} function.
#' Its Help file should be consulted for its full argument list.
#'
#' @param abbrevLength Default = 0. A non-zero value of this argument is used as
#' the \code{minlength} argument of the \code{\link{abbreviate}} function to
#' abbreviate the names of the conditioning factors both in the automatically
#' generated console legend and the optional plot legend.
#'
#' @param ... Further arguments to pass down to either the \code{legend} function
#' or \code{\link[lattice]{print.trellis}}, which should be consulted for details.
#' Care should be taken to ensure that the names of arguments do not coincide.
#' See \code{\link{defaultStrucLegend}} for arguments for the default \code{legend}
#' function.
#'
#' @note Do not use the \code{packet.panel} argument of \code{print.trellis},
#'  as this will totally mess up the display.
#'
#' @seealso \code{\link[lattice]{print.trellis}}
#'
#' @examples
#' require(datasets)
#' # quakes data
#' #
#' # Create and save plot
#' out <- strucplot(lat ~ long|cut(mag,5)*cut(depth,4), data = quakes,
#'   col="blue", main = "Earthquake locations, by magnitude and depth")
#'
#' # Summary:
#' summary(out)
#'
#' # Default output -- structure legend on console only
#'    print(out)
#'
#' # Add legends to the plot on either right or bottom (note partial matching)
#'    print(out, legendLoc = "right")
#'    print(out, legendLoc = "b")
#' #
#' # Abbreviate the conditioning factor names
#'    print(out, legendLoc = "b", abbrev = 5)
#' #
#' # Plot the legend by itself on a separate page
#'    print(out, legendLoc = "newp")
#' #
#' # Extra grid "gp" arguments to alter text appearance
#'    print(out, legendLoc = "b",col="blue",fontface = "italic",
#'    abbrev = 5)
#'
print.structured <- function(x
  ,legendLoc = c("left","right","top","bottom","newpage")
  ## optional legend location for a plotted legend
  ,legend = defaultStrucLegend
  ,abbrevLength = 0 ## minlength argument for abbreviating conditioning variable
#       names
   ,...  ## further named arguments to legend() or print.trellis()
  )
{
  struc <- attr(x,"structure")
## create text string for console display
  txt <-displayStruc(struc= struc ,abbrevLength = abbrevLength)
  ## and display
  cat(sprintf("\nPLOT STRUCTURE\n\n%s\n\n%s",txt[1],txt[2]))
  ## create legend if needed
  if(!missing(legendLoc)){
    pos <- tryCatch(match.arg(legendLoc),
                    error = function(e)stop("Incorrect legendLoc argument"))
    leg <- legend(attr(x,"structure"),
                  legendLoc = pos,
                  abbrevLength = abbrevLength,
                  ...)
    if(!pos == "newpage"){
      oleg <- x$legend
      leg <- structure(list(list(fun = leg)), names = pos)
      if(pos %in% names(oleg))
        stop("Structure legend position conflicts with plot legend position")
        else x$legend <- c(oleg,leg)
    }
    NextMethod()
    if(pos == "newpage"){
        ans <- readline("\nHit <return> to show structure legend:")
        grid.newpage()
        grid.draw(leg)
      }
  } else NextMethod()
}

#
#' @rdname print.structured
plot.structured <- function(...) print(...)
#
#' @rdname  print.structured
summary.structured <- function(object,...)
{
  out <- list(
  call = object[["call"]],
  formula = attr(object,"formula"),
  structure = attr(object,"structure"))
  class(out) <- c("summary.structured", "list")
  out
}
#
#' @rdname  print.structured
print.summary.structured <- function(x,...)
{
  cat("\n  xyplot call:\n")
  print(x$call)
  cat("\n Actual formula:\n")
  print(x$formula)
  displayStruc(x$structure)
}

