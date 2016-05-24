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

"xpose.bin" <- function (data,
                         y,
                         bins=10) {

    ## substitute for equal count algorithm in trellis 
    ## bin a continuous variable for a bwplot
    if (length(unique(data[[y]])) >= bins) {
    for (i in 1:length(names(data))) {
      if (names(data)[i] == y) {
        y.data <- data[i]
        mxr <- max(data[i])
        mnr <- min(data[i])
      }
      #if (names(data)[i] == x) {
      #  x.data <- data[i]
      #}
    }

    mdif <- mxr - mnr
    mit  <- mdif/(bins)
    curr <- mnr
    
    
    binlist <- c()
    for (i in 1:bins) {
      binlist <- c(binlist, curr)
      curr <- curr + mit
      if (i == bins) {
        binlist <- c(binlist, curr)
      }
    }
    #colnames(x.data) <- c("x")
    colnames(y.data) <- c("y")
    ## now assign all the bits of data
    bwdata <- cut(y.data$y, binlist, include.lowest = T)

    } else {
      # categorical
      for (i in 1:length(names(data))) {
        if (names(data)[i] == y) {
          y.data <- data[i]
        }
      }
      binlist <- c()
      for (i in 1:length(unique(y))) {
        binlist <- c(binlist, y[i])
        if (i == bins) {
          binlist <- c(binlist, y[i])
        }
      }
    #colnames(x.data) <- c("x")
    colnames(y.data) <- c("y")
    ## now assign all the bits of data
    bwdata <- y.data$y
    }
    

    return(bwdata)
}