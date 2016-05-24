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


#' @title Read the Turtle's Status
#'
#' @description
#' This function gives information about the current
#' Turtle's position, direction, and on display options.
#' 
#' @details
#' The Turtle must be initialized prior to using
#' this function, see \code{\link{turtle_init}}. 
#' 
#' @return
#' Returns a list with three elements.
#' 
#' @family TurtleGraphics
#' @rdname turtle_status
#' @export
turtle_status <- function()
{
   .turtle_check()
   
   gp <- get("gpar_path", envir=.turtle_data)
   
   list(DisplayOptions=list(
         col=gp$col,
         lty=gp$lty,
         lwd=gp$lwd,
         visible=get("visible", envir=.turtle_data),
         draw=get("draw", envir=.turtle_data)
      ),
      Terrarium=list(
         width=get("width", envir=.turtle_data),
         height=get("height", envir=.turtle_data)
      ),
      TurtleStatus=list(
         x=get("x", envir=.turtle_data),
         y=get("y", envir=.turtle_data),
         angle=get("angle", envir=.turtle_data)
      ))
}
