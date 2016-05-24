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


# This object shall not be exported:
.turtle_data <- new.env()


# This function shall not be exported:
.turtle_check <- function()
{
   if (!exists("width", envir=.turtle_data))
      stop("The turtle has not been initialized, please call turtle_init() first.")
}


# This function shall not be exported:
.turtle_undraw <- function()
{
   grid.remove("turtle_head")
   grid.remove("turtle_body")
   grid.remove("turtle_leg1")
   grid.remove("turtle_leg2")
   grid.remove("turtle_leg3")
   grid.remove("turtle_leg4")
   invisible(NULL)
}


# This function shall not be exported:
.turtle_draw <- function()
{
   # current turtle status:
   x <- get("x", envir=.turtle_data)
   y <- get("y", envir=.turtle_data)
   angle <- pi * get("angle", envir=.turtle_data) / 180.0
   
   gp1 <- get("gpar_turtle1", envir=.turtle_data)
   gp2 <- get("gpar_turtle2", envir=.turtle_data)
   
   scale <- max(get("width", envir=.turtle_data), get("height", envir=.turtle_data))
   tmp <- 0.03*scale
   grid.circle(
      x + tmp * sin(angle),
      y + tmp * cos(angle),
      default.units='native',  
      gp = gp1,
      r=0.015*scale,
      name = "turtle_head")
   grid.circle(
      x + tmp * sin(angle + pi/3),
      y + tmp * cos(angle + pi/3),
      default.units='native',  
      gp = gp1,
      r=0.01*scale,
      name = "turtle_leg1")
   grid.circle(
      x + tmp * sin(angle - pi/3),
      y + tmp * cos(angle - pi/3),
      default.units='native',  
      gp = gp1,
      r=0.01*scale, name = "turtle_leg2")
   grid.circle(
      x + tmp * sin(angle + 2*pi/3),
      y + tmp * cos(angle + 2*pi/3),  
      gp = gp1,
      r=0.01*scale, 
      default.units='native',
      name = "turtle_leg3")
   grid.circle(
      x + tmp * sin(angle - 2*pi/3),
      y + tmp * cos(angle - 2*pi/3),
      default.units='native',  
      gp = gp1,
      r=0.01*scale,
      name = "turtle_leg4")
   grid.circle(
      x,
      y,
      r=0.03*scale,
      default.units='native',
      gp = gp2,
      name = "turtle_body")
   invisible(NULL)
}


# This function shall not be exported:
.turtle_redraw <- function()
{
   if (get("visible", envir=.turtle_data)) {
      .turtle_undraw()
      .turtle_draw()
   }
   invisible(NULL)
}


# This function shall not be exported:
.turtle_set_default_params <- function()
{
   # changes neither height nor width,
   # which is done by turtle_init()
   assign(envir=.turtle_data, "gpar_path",
      gpar(col="black", lty=1, lwd=1))
   assign(envir=.turtle_data, "gpar_turtle1",
      gpar(fill = "darkorange2", col = "brown"))
   assign(envir=.turtle_data, "gpar_turtle2",
      gpar(fill = "darkorange4", col = "brown"))
   assign(envir=.turtle_data, "visible",  TRUE)
   assign(envir=.turtle_data, "draw",     TRUE)
   assign(envir=.turtle_data, "x",        get("width", envir=.turtle_data)*0.5)
   assign(envir=.turtle_data, "y",        get("height", envir=.turtle_data)*0.5)
   assign(envir=.turtle_data, "angle",    0.0)   
}
