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


#' @title Move the Turtle Forward or Backward
#'
#' @description
#' \code{turtle_forward} moves the Turtle in forward direction and
#' \code{turtle_backward} moves the Turtle back.
#' 
#' 
#' @details
#' The Turtle must be initialized prior to using
#' these functions, see \code{\link{turtle_init}}. 
#' 
#' These functions make use of the Turtle's display options specified by
#' the \code{\link{turtle_param}} function
#' (or if not, use the default options set by \code{\link{turtle_init}}).
#' 
#' Note that if \code{\link{turtle_up}} or \code{\link{turtle_down}}
#' was called, the Turtle's trace will be or not be drawn, respectively.
#' 
#' If you are willing to call these functions in an R loop,
#' you may want to hide the Turtle temporarily (see \code{\link{turtle_hide}}
#' and \code{\link{turtle_do}})
#' before making actual moves. This will increase the drawing performance
#' significantly.
#' 
#' 
#' @param distance single numeric value; specifies the distance to make.
#'        Negative distance results in moving in the opposite direction.
#' @param direction character string; moving direction.
#'        One of \code{"forward"} or \code{"backward"}.
#' 
#'
#'
#' @examples
#' turtle_init()
#' turtle_left(30)
#' turtle_forward(2)
#' turtle_up()
#' turtle_forward(1)
#' turtle_down()
#' turtle_right(60)
#' turtle_forward(9)
#' 
#' @family TurtleGraphics
#' @rdname turtle_move
#' @export
turtle_move <- function(distance, direction = c("forward", "backward"))
{
   .turtle_check()
   
   direction <- match.arg(direction)
   stopifnot(is.numeric(distance), length(distance) == 1, is.finite(distance))
   
   if (distance < 0)
      warning("Negative value of `distance` moves the Turtle in the opposite direction.")
  
   if (direction == 'backward')
      distance <- -distance
   
   # current values for .turtle_history
   curX        <- get("x", envir=.turtle_data)
   curY        <- get("y", envir=.turtle_data)
   curAng      <- get("angle", envir=.turtle_data)
   curGp       <- get("gpar_path", envir=.turtle_data)
   curDraw     <- get("draw", envir=.turtle_data)
   curMode     <- get("mode", envir=.turtle_data)
   curWidth    <- get("width", envir=.turtle_data)
   curHeight   <- get("height", envir=.turtle_data)
   curVisible  <- get("visible", envir=.turtle_data)
   
   if (curMode == 'error')
      .turtle_draw_error(distance, curX, curY, curAng, curGp, curDraw,
         curWidth, curHeight, curVisible)
   else if (curMode == 'clip')
      .turtle_draw_clip(distance, curX, curY, curAng, curGp, curDraw,
         curWidth, curHeight, curVisible)
   else
      .turtle_draw_cycle(distance, curX, curY, curAng, curGp, curDraw,
         curWidth, curHeight, curVisible)
   
   invisible(NULL)
}


#' @rdname turtle_move
#' @export
turtle_forward <- function(distance)
{
   turtle_move(distance, 'forward')
}


#' @rdname turtle_move
#' @export
turtle_backward <- function(distance)
{
   turtle_move(distance, 'backward')
}


# This function shall not be exported:
.turtle_draw_fromto <- function(curX, curY, newX, newY, curGp, curDraw, curVisible)
{
   if (curVisible) .turtle_undraw()
   
   if (curDraw) {
      grid.lines(c(curX, newX), c(curY, newY), 
         gp = curGp, default.units='native')
   }
   
   assign("x", newX, envir=.turtle_data)
   assign("y", newY, envir=.turtle_data)
   
   if (curVisible) .turtle_draw()
}


# This function shall not be exported:
.turtle_draw_error <- function(distance, curX, curY, curAng, curGp, curDraw,
         curWidth, curHeight, curVisible)
{
   newX <- curX + distance * sin(curAng * pi / 180)
   newY <- curY + distance * cos(curAng * pi / 180)
   
   if (newY > curHeight || newY < 0 || newX > curWidth || newX < 0)
      stop("The Turtle escaped from the terrarium. :-(")
   
   .turtle_draw_fromto(curX, curY, newX, newY, curGp, curDraw, curVisible)
}


# This function shall not be exported:
.turtle_draw_clip <- function(distance, curX, curY, curAng, curGp, curDraw,
         curWidth, curHeight, curVisible)
{
   # new position
   newX <- curX + distance * sin(curAng * pi / 180)
   newY <- curY + distance * cos(curAng * pi / 180)
   
   .turtle_draw_fromto(curX, curY, newX, newY, curGp, curDraw, curVisible)
}


# This function shall not be exported:
.turtle_draw_cycle <- function(distance, curX, curY, curAng, curGp, curDraw,
                               curWidth, curHeight, curVisible)
{
  tmpVisible <- curVisible
  if (curVisible){
    tmpVisible <- TRUE
    curVisible <- FALSE
    .turtle_undraw()
    assign("visible", FALSE, envir=.turtle_data)
  }  
  newX <- curX + distance * sin(curAng * pi / 180)
  newY <- curY + distance * cos(curAng * pi / 180)  
  tmpAng <- curAng %% 360
  
  # if turtle leaves terrarium:
  if(newX > curWidth | newX < 0 | newY > curHeight | newY < 0){
    
    if(tmpAng < 90){
      distWallX <- curWidth - curX
      distWallY <- curHeight - curY
      if(distWallX / distWallY > tan(tmpAng * pi / 180)){
        # up
        .turtle_cycle_up(distance, curX, curY, curAng, curGp, curDraw,
                         curWidth, curHeight, curVisible,
                         newX, newY, tmpAng)
      }else{
        # right
        .turtle_cycle_right(distance, curX, curY, curAng, curGp, curDraw,
                            curWidth, curHeight, curVisible,
                            newX, newY, tmpAng)      
      }
    }else{
      if(tmpAng < 180){
        distWallX <- curWidth - curX
        distWallY <- curY
        if(distWallY / distWallX > tan((tmpAng - 90) * pi / 180)){
          # right
          .turtle_cycle_right(distance, curX, curY, curAng, curGp, curDraw,
                              curWidth, curHeight, curVisible,
                              newX, newY, tmpAng)
        }else{
          # down
          .turtle_cycle_down(distance, curX, curY, curAng, curGp, curDraw,
                             curWidth, curHeight, curVisible,
                             newX, newY, tmpAng)      
        }
        
      }else{
        if(tmpAng < 270){
          distWallX <- curX
          distWallY <- curY
          if(distWallX / distWallY > tan((tmpAng - 180) * pi / 180)){
            # down
            .turtle_cycle_down(distance, curX, curY, curAng, curGp, curDraw,
                               curWidth, curHeight, curVisible,
                               newX, newY, tmpAng)
          }else{
            # left
            .turtle_cycle_left(distance, curX, curY, curAng, curGp, curDraw,
                               curWidth, curHeight, curVisible,
                               newX, newY, tmpAng)      
          }
          
        }else{
          distWallX <- curX
          distWallY <- curHeight - curY
          if(distWallY / distWallX > tan((tmpAng - 270) * pi / 180)){
            # left
            .turtle_cycle_left(distance, curX, curY, curAng, curGp, curDraw,
                               curWidth, curHeight, curVisible,
                               newX, newY, tmpAng)
          }else{
            # up
            .turtle_cycle_up(distance, curX, curY, curAng, curGp, curDraw,
                             curWidth, curHeight, curVisible,
                             newX, newY, tmpAng)      
          }
        }
      }      
    }
  }
  else .turtle_draw_fromto(curX, curY, newX, newY, curGp, curDraw, curVisible)

  if (tmpVisible){
    .turtle_draw()
    assign("visible", TRUE, envir=.turtle_data)
  }
}


# This function shall not be exported:
.turtle_cycle_up <- function(distance, curX, curY, curAng, curGp, curDraw,
                             curWidth, curHeight, curVisible,
                             newX, newY, tmpAng)
{
  distPassed <- (curHeight - curY) / cos(tmpAng * pi / 180)
  nextX <- curX + distPassed * sin(tmpAng * pi / 180)
  distLeft <- distance - distPassed
  
  .turtle_draw_fromto(curX, curY, nextX, curHeight, curGp, curDraw, curVisible)
  assign("y", 0, envir=.turtle_data)      
  .turtle_draw_cycle(distLeft, nextX, 0, curAng, curGp, curDraw,
                     curWidth, curHeight, curVisible)
}


# This function shall not be exported:
.turtle_cycle_down <- function(distance, curX, curY, curAng, curGp, curDraw,
                               curWidth, curHeight, curVisible,
                               newX, newY, tmpAng)
{
  tmpAng <- tmpAng - 180
  
  distPassed <- (curY) / cos(tmpAng * pi / 180)
  nextX <- curX - distPassed * sin(tmpAng * pi / 180)
  distLeft <- distance - distPassed
  
  .turtle_draw_fromto(curX, curY, nextX, 0, curGp, curDraw, curVisible)
  assign("y", curHeight, envir=.turtle_data)      
  .turtle_draw_cycle(distLeft, nextX, curHeight, curAng, curGp, curDraw,
                     curWidth, curHeight, curVisible)
}


# This function shall not be exported:
.turtle_cycle_left <- function(distance, curX, curY, curAng, curGp, curDraw,
                               curWidth, curHeight, curVisible,
                               newX, newY, tmpAng)
{
  tmpAng <- tmpAng - 270
  
  distPassed <- (curX) / cos(tmpAng * pi / 180)
  nextY <- curY + distPassed * sin(tmpAng * pi / 180)
  distLeft <- distance - distPassed
  
  .turtle_draw_fromto(curX, curY, 0, nextY, curGp, curDraw, curVisible)
  assign("x", curWidth, envir=.turtle_data)
  .turtle_draw_cycle(distLeft, curWidth, nextY, curAng, curGp, curDraw,
                     curWidth, curHeight, curVisible)  
}


# This function shall not be exported:
.turtle_cycle_right <- function(distance, curX, curY, curAng, curGp, curDraw,
                                curWidth, curHeight, curVisible,
                                newX, newY, tmpAng)
{
  tmpAng <- tmpAng - 90

  distPassed <- (curWidth - curX) / cos(tmpAng * pi / 180)
  nextY <- curY - distPassed * sin(tmpAng * pi / 180)
  distLeft <- distance - distPassed
  
  .turtle_draw_fromto(curX, curY, curWidth, nextY, curGp, curDraw, curVisible)
  assign("x", 0, envir=.turtle_data)
  .turtle_draw_cycle(distLeft, 0, nextY, curAng, curGp, curDraw,
                     curWidth, curHeight, curVisible) 
}
