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


#' @title
#' Set the Turtle's Position and Direction
#'
#' @description
#' \code{turtle_goto} and \code{turtle_setpos} move the Turtle to a 
#' given location without changing its direction.
#' 
#' \code{turtle_setangle} 
#' rotates the Turtle to a given (absolute) angle,
#' where 0 denotes a north-facing Turtle.
#'
#' @details
#' The terrarium must be initialized prior to using
#' these functions, see \code{\link{turtle_init}}.
#' 
#' If the given location (x, y) lies outside the terrarium,
#' the behavior of these functions depends
#' on the \code{mode} argument in \code{\link{turtle_init}}.
#' 
#' \code{turtle_goto} may draw the path between the current Turtle's 
#' position and the new location. Its behavior depends on the current
#' plot settings, 
#' see \code{\link{turtle_up}}, \code{\link{turtle_down}}. In case of 
#' \code{turtle_setpos}, however, the path drawing is always disabled.
#'   
#' @param x,y  numeric; coordinates specifying new Turtle's location.
#' @param angle numeric; rotation angle in degrees.
#'
#' @family TurtleGraphics
#' @rdname turtle_goto
#' @export
turtle_goto <- function(x, y)
{
   stopifnot(is.numeric(x), length(x)==1, is.finite(x))
   stopifnot(is.numeric(y), length(y)==1, is.finite(y))

   .turtle_check()
   
   curX        <- get("x", envir=.turtle_data)
   curY        <- get("y", envir=.turtle_data)
   curAng      <- get("angle", envir=.turtle_data)
   curDraw     <- get("draw", envir=.turtle_data)
   curVisible  <- get("visible", envir=.turtle_data)
   curGp       <- get("gpar_path", envir=.turtle_data)
   curMode     <- get("mode", envir=.turtle_data)
   curWidth    <- get("width", envir=.turtle_data)
   curHeight   <- get("height", envir=.turtle_data)
   
   if (curMode == 'clip' ||
     (x >= 0 && x<= curWidth && y>=0 && y<=curHeight)) {
      .turtle_draw_fromto(curX, curY, x, y, curGp, curDraw, curVisible)      
   }
   else {
      stop('Given coordinates lie outside the terrarium. :-(')
   }
   
   invisible(NULL)
}


#' @rdname turtle_goto
#' @export
turtle_setpos <- function(x, y)
{
   stopifnot(is.numeric(x), length(x)==1, is.finite(x))
   stopifnot(is.numeric(y), length(y)==1, is.finite(y))
   
   .turtle_check()
   
   curMode     <- get("mode", envir=.turtle_data)
   curWidth    <- get("width", envir=.turtle_data)
   curHeight   <- get("height", envir=.turtle_data)
   
   
   if (curMode == 'clip' || (x >= 0 && x<= curWidth && y>=0 && y<=curHeight)) {
      assign("x", x, envir=.turtle_data)
      assign("y", y, envir=.turtle_data)
      .turtle_redraw()      
   }
   else {
      stop('Given coordinates lie outside the terrarium. :-(')
   }
   
   invisible(NULL)
}


#' @rdname turtle_goto
#' @export
turtle_setangle <- function(angle)
{
   stopifnot(is.numeric(angle), length(angle)==1, is.finite(angle))
   
   .turtle_check()
   assign("angle", angle, envir=.turtle_data)
   .turtle_redraw()      

   invisible(NULL)
}
