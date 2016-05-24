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
#' Turn (Rotate) the Turtle
#'
#' @description
#' Turn the Turtle in the given direction by the given angle.
#' 
#' @details
#' The Turtle must be initialized prior to using
#' this function, see \code{\link{turtle_init}}. 
#' 
#' 
#' @param angle  single numeric value; rotation angle in degrees.
#'    A negative value turns the Turtle in the opposite
#'    direction than the given one.
#' @param direction character string; direction of the turn.
#'    Possible values are \code{"left"} and \code{"right"}.
#'
#' @examples
#' turtle_init()
#' turtle_left(30) # equivalent to turtle_turn(30, "left")
#' turtle_right(40)
#' turtle_turn(30, sample(c("left", "right"), 1)) # random turn
#' 
#' @family TurtleGraphics
#' @export
#' @rdname turtle_turn
turtle_turn <- function(angle, direction = c("left", "right"))
{
   .turtle_check()
   
   direction <- match.arg(direction)
   stopifnot(is.numeric(angle), length(angle) == 1, is.finite(angle))
   
   if (angle < 0)
      warning("Negative `angle` turns the Turtle in the opposite direction.")
   
   if (direction == "left")
      angle = -angle
   
   assign(envir=.turtle_data, "angle", get("angle", envir=.turtle_data) + angle)
   
   .turtle_redraw()
   
   invisible(NULL)
}


#' @rdname turtle_turn
#' @export
turtle_left <- function(angle)
{
   turtle_turn(angle = angle, direction = "left")
}


#' @rdname turtle_turn
#' @export
turtle_right <- function(angle)
{
   turtle_turn(angle = angle, direction = "right")
}
