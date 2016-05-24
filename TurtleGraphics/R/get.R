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
#' Get the Turtle's Current Position and Direction
#'
#' @description
#' \code{turtle_getpos} returns the Turtle's current position on the plane.
#' 
#' \code{turtle_getangle} returns the Turtle's current direction,
#' in degrees. An angle of 0 represents a north-facing Turtle.
#' 
#' @details
#' The terrarium must be initialized prior to using these functions,
#' see \code{\link{turtle_init}}.
#' 
#' 
#' 
#' @return
#' Both functions return a (named) numeric vector.
#' \code{turtle_getpos} returns a vector of length two which
#' specifies the \code{x} and \code{y} coordinates.
#' The \code{turtle_getangle} returns the angle.
#'    
#' @examples
#' turtle_init()
#' turtle_getpos()["x"] # x coordinate
#' turtle_getpos()["y"] # y coordinate
#'  
#' @family TurtleGraphics
#' @rdname turtle_getpos
#' @export
turtle_getpos <- function()
{
   .turtle_check()
   c(x=get("x", envir=.turtle_data),
     y=get("y", envir=.turtle_data))
}

#' @rdname turtle_getpos
#' @export
turtle_getangle <- function()
{
   .turtle_check()
   c(angle=get("angle", envir=.turtle_data))
}
