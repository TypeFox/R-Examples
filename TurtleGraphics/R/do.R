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
#' Evaluate a Larger Portion of Turtle Drawing Code
#'
#' @description
#' \code{turtle_do} evaluates an R expression 
#' with the Turtle temporarily hidden (for performance reasons).
#' 
#' @param expr expression to evaluate
#' 
#'  
#' @details
#' The terrarium must be initialized prior to using
#' these functions, see \code{\link{turtle_init}}.
#' 
#' In order to decrease the evaluation time of \code{expr},
#' it is evaluated with Turtle temporarily hidden.
#' Basically it means that if a Turtle image is visible (see 
#' \code{\link{turtle_show}} and \code{\link{turtle_hide}}) \code{turtle_do}
#' removes it, evaluates \code{expr} and redraws it on the function exit.
#'    
#' @examples
#' turtle_init()
#' turtle_do({
#'    for (i in 1:4) {
#'       turtle_forward(50)
#'       turtle_right(90)
#'    }
#' })
#' 
#' @family TurtleGraphics
#' @rdname turtle_do
#' @export
turtle_do <- function(expr){
   .turtle_check()
   curVisible <- get("visible", envir=.turtle_data)
   if (curVisible)
      turtle_hide()
   
   eval(substitute(expr), enclos = parent.frame())
   
   if (curVisible)
      turtle_show()
 
   invisible(NULL)
}
   