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


#' @title Show or Hide the Turtle
#'
#' @description
#' These functions enable or disable displaying the Turtle's
#' image on the screen.
#' 
#' @details
#' The Turtle must be initialized prior to using
#' this function, see \code{\link{turtle_init}}. 
#' 
#' It is recommended to hide the Turtle when
#' performing multiple Turtle moves, for efficiency reasons,
#' see also \code{\link{turtle_do}}.
#'
#' @examples
#' turtle_init()
#' turtle_forward(4)
#' turtle_hide()
#' turtle_left(30)
#' turtle_forward(3)
#' 
#' @family TurtleGraphics
#' @rdname turtle_show
#' @export
turtle_show <- function()
{
   .turtle_check()
   
   if (!get("visible", envir=.turtle_data)) {
      assign(envir=.turtle_data, "visible", TRUE)
      .turtle_draw()      
   }
   else{
      warning("The turtle is already visible.")
   }
   
   invisible(NULL)
}


#' @rdname turtle_show
#' @export
turtle_hide <- function()
{
   .turtle_check()
   
   if (get("visible", envir=.turtle_data)) {
      assign(envir=.turtle_data, "visible", FALSE)
      .turtle_undraw()
   }
   else{
      warning("The turtle is already hidden.")
   }
   
   invisible(NULL)
}
