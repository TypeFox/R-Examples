# Xpose 4
# An R-based population pharmacokinetic/
# pharmacodynamic model building aid for NONMEM.
# Copyright (C) 1998-2004 E. Niclas Jonsson and Mats Karlsson.
# Copyright (C) 2005-2008 Andrew C. Hooker, Justin J. Wilkins, 
# Mats O. Karlsson and E. Niclas Jonsson.
# Copyright (C) 2009-2010 Andrew C. Hooker, Mats O. Karlsson and 
# E. Niclas Jonsson.

# This file is a part of Xpose 4.
# Xpose 4 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation, either version 3
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with this program.  A copy can be cound in the R installation
# directory under \share\licenses. If not, see http://www.gnu.org/licenses/.

"fix.wrapped.lines" <- function(par.mat)
{
  ## This function unwraps matrix lines that NONMEM and nmsee wraps. Assumes
  ## no more than 60 etas
  nlines <- length(par.mat)
  unw.par <- list()
  unw.par[[1]] <- par.mat[[1]]
  ## first line is never wrapped!
  if(length(par.mat) == 1) return(unw.par)
  unw.line <- 2
  doneline <- 0
  doneline1 <- 0
  doneline2 <- 0
  doneline3 <- 0
  for(i in 2:nlines) {
    if(unw.line <= 12) {
      unw.par[[unw.line]] <- par.mat[[i]]
      unw.line <- unw.line + 1
    } else if(unw.line > 12 && unw.line <= 24) {
      if(doneline == i) {
        doneline <- 0
        next
      }
      unw.par[[unw.line]] <- c(par.mat[[i]], par.mat[[i + 1]])
      unw.line <- unw.line + 1
      doneline <- i + 1
    }else if(unw.line > 24 && unw.line <= 36) {
      if(doneline == i) {
        doneline <- 0
        next
      }
      if(doneline1 == i) {
        doneline1 <- 0
        next
      }
      unw.par[[unw.line]] <- c(par.mat[[i]],
                               par.mat[[i +1]],
                               par.mat[[i + 2]])
      unw.line <- unw.line + 1
      doneline <- i + 1
      doneline1 <- i + 2
    }else if(unw.line > 36 && unw.line <= 48) {
      if(doneline == i) {
        doneline <- 0
        next
      }
      if(doneline1 == i) {
        doneline1 <- 0
        next
      }
      if(doneline2 == i) {
        doneline2 <- 0
        next
      }
      unw.par[[unw.line]] <- c(par.mat[[i]],
                               par.mat[[i +1]],
                               par.mat[[i +2]],
                               par.mat[[i +3]])
      unw.line <- unw.line + 1
      doneline <- i + 1
      doneline1 <- i + 2
      doneline2 <- i + 3
    } else if(unw.line > 48 && unw.line <= 60) {
      if(doneline == i) {
        doneline <- 0
        next
      }
      if(doneline1 == i) {
        doneline1 <- 0
        next
      }
      if(doneline2 == i) {
        doneline2 <- 0
        next
      }
      if(doneline3 == i) {
        doneline3 <- 0
        next
      }
      unw.par[[unw.line]] <- c(par.mat[[i]],
                               par.mat[[i +1]],
                               par.mat[[i +2]],
                               par.mat[[i +3]],
                               par.mat[[i +4]])
      unw.line <- unw.line + 1
      doneline <- i + 1
      doneline1 <- i + 2
      doneline2 <- i + 3
      doneline3 <- i + 4
    }
  }
  return(unw.par)
}
