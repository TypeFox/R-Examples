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


#' @title Reset the Turtle's Position and Direction
#'
#' @description
#' This function resets the Turtle's position, direction,
#' and graphical options.
#' 
#' @details
#' The Turtle must be initialized prior to using
#' this function, see \code{\link{turtle_init}}. 
#' 
#' After a call to this function, the Turtle will be placed
#' in the terrarium's center and it will be directed
#' to the north.
#'
#' The drawing remains unchanged.
#' 
#' 
#' @examples
#' turtle_init()
#' turtle_forward(4)
#' turtle_param(col="red", lty=2, lwd=3)
#' turtle_reset()
#' turtle_left(45)
#' turtle_forward(3)
#' 
#' @family TurtleGraphics
#' @export
#' @rdname turtle_reset
turtle_reset <- function()
{
   .turtle_check()
   
   was_visible <- get("visible", envir=.turtle_data)
   
   .turtle_set_default_params()
  
   if (!was_visible)
      .turtle_draw()
   else
      .turtle_redraw()
   
   invisible(NULL)
} 
