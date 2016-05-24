# 
#  Algorithm for computing exact Oja Median by bounded search.
#  Copyright (C) 2015 Oleksii Pokotylo
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful, 
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
`ojaMedianExB` <- function(X, control=ojaMedianControl(...), ...){
  
  x <- y <- 1
  rows <- dim(X)[1]
  cols <- dim(X)[2]
  outvec <- c(1:cols)
  storage.mode(rows) <- "integer"
  storage.mode(cols) <- "integer"
  storage.mode(X) <- "double"
  storage.mode(outvec) <- "double"
  
  action <- 6
  param2 <- control$volume
  param3 <- control$boundedExact  # if FALSE the approximative method is used
  param4 <- debug <- 0
  # debug = 1

  res<-.C("r_oja", rows, cols, X, vec = outvec, y, as.integer(action), as.double(control$maxlines), as.double(param2), as.integer(param3), as.integer(param4), as.integer(debug),1)
  RES <- res$vec

  names(RES)<-colnames(X)
  return(RES)
}
