{###############################################################################
# consecutive.R
# This file is part of the R package dostats.
# 
# Copyright 2012 Andrew Redd
# Date: 6/1/2012
# 
# DESCRIPTION
# ===========
# computes a vector that changes every time the element is different from the 
# previous.
# 
# LICENSE
# ========
# dostats is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software 
# Foundation, either version 3 of the License, or (at your option) any later 
# version.
# 
# dostats is distributed in the hope that it will be useful, but WITHOUT ANY 
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with 
# dostats. If not, see http://www.gnu.org/licenses/.
# 
}###############################################################################

#' compute an indicator to group consecutive values
#'
#' computes a vector that changes every time the element is different from the 
#' previous.
#'
#' @param x a vector
#' @param ... ignored, might be used for forward compatibility.
#' @return an integer vector.
#' @export
seq_consecutive <-
function(x,...){
  stopifnot(is.vector(x))
  y <- integer(length(x))
  y[1] <- 1
  if(length(x)>1)for(i in seq_len(length(x))[-1]){
    if(x[i]==x[i-1])
      y[i] <- y[i-1]
    else
      y[i] <- y[i-1] +1
  }
  return(y)
}
