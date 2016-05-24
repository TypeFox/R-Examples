##    TurtleGraphics package for R
##    Copyright (C) 2014 Rexamine
##
##    This program is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program.  If not, see <http://www.gnu.org/licenses/>.


#' @title Set Display Options
#'
#' @description
#' Sets the display options for the Turtle's trace.
#' It is possible to change its color, line type and
#' line width.
#' 
#'   
#' @param col  numeric or character; trace color, see e.g. \code{\link{colors}}
#'      and \code{\link{gpar}}.
#' @param lty  numeric; trace line type, see \code{\link{gpar}}.
#' @param lwd  numeric; trace line width, see \code{\link{gpar}}.
#' 
#' @details
#' The Turtle must be initialized prior to using
#' this function, see \code{\link{turtle_init}}. 
#' 
#' @examples
#' turtle_init()
#' turtle_forward(5)
#' turtle_up()
#' turtle_forward(3)
#' turtle_down()
#' turtle_left(90)
#' turtle_forward(5)
#' turtle_param(col = "red", lwd = 2, lty = 2)
#' turtle_forward(5)
#' 
#' @family TurtleGraphics
#' @rdname turtle_param
#' @export
turtle_param <- function(col=NULL, lwd=NULL, lty=NULL)
{
   .turtle_check()
   
   if (all(is.null(col), is.null(lwd), is.null(lty)))
      stop("You need to provide at least one of: `col`, `lwd`, `lty`.")
   
   if (!is.null(col)) {
      tryCatch({
         tmp <- col2rgb(col)
      },
      error = function(e) {
         stop("Given `col` is not a correct color specifier.")
      })
      
      gp <- get("gpar_path", envir=.turtle_data)
      gp$col <- col
      assign(envir=.turtle_data, "gpar_path", gp)
   }
   
   if (!is.null(lty)) {
      stopifnot(is.numeric(lty), length(lty) == 1, is.finite(lty), lty > 0)
      gp <- get("gpar_path", envir=.turtle_data)
      gp$lty <- lty
      assign(envir=.turtle_data, "gpar_path", gp)
   }
   
   if (!is.null(lwd)) {
      stopifnot(is.numeric(lwd), length(lwd) == 1, is.finite(lwd), lwd > 0)
      gp <- get("gpar_path", envir=.turtle_data)
      gp$lwd <- lwd
      assign(envir=.turtle_data, "gpar_path", gp)
   }
  
   invisible(NULL)
}


#' @rdname turtle_param
#' @export
turtle_col <- function(col)
{
   turtle_param(col = col)
}


#' @rdname turtle_param
#' @export
turtle_lwd <- function(lwd)
{
   turtle_param(lwd = lwd)
}


#' @rdname turtle_param
#' @export
turtle_lty <- function(lty)
{
   turtle_param(lty = lty)
}



#' @title
#' Turn on or off Turtle Trace Drawing
#' 
#' @description
#' When the Turtle moves, it may or may not
#' leave a visible trace. These functions
#' control such a behavior.
#' 
#' @details
#' The Turtle must be initialized prior to using
#' this function, see \code{\link{turtle_init}}. 
#' 
#' @family TurtleGraphics
#' @rdname turtle_up
#' @export
turtle_up <- function()
{
   .turtle_check()
   assign(envir=.turtle_data, "draw", FALSE)
   invisible(NULL)
}


#' @rdname turtle_up
#' @export
turtle_down <- function()
{
   .turtle_check()
   assign(envir=.turtle_data, "draw", TRUE)
   invisible(NULL)
}
