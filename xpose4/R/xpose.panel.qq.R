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

xpose.panel.qq <- function(x, object,
                                  #subscripts,
                                  #inclZeroWRES = FALSE,
                                  #onlyfirst = FALSE,
                                  #samp = NULL,
                                  #xvarnam=NULL,
                                  #breaks=NULL,
                                  #ylab=ylb,
                                  #xlab=xlb,
                                  pch=object@Prefs@Graph.prefs$pch,
                                  col=object@Prefs@Graph.prefs$col,
                                  cex=object@Prefs@Graph.prefs$cex,
                                  abllty = object@Prefs@Graph.prefs$abllty,
                                  abllwd = object@Prefs@Graph.prefs$abllwd,
                                  ablcol = object@Prefs@Graph.prefs$ablcol,
                                  grid = object@Prefs@Graph.prefs$grid,
                                  ...) {

#  browser()
##   if(!is.null(samp)) {
##     data <- SData(object,inclZeroWRES,onlyfirst=onlyfirst,samp=samp)
##   } else {
##     data <- Data(object,inclZeroWRES,onlyfirst=onlyfirst)
##   }
##   if(length(unique(data[subscripts,xvarnam])) <= object@Prefs@Cat.levels)
##     x <- as.factor(x)
    
  if(grid != FALSE) {
    panel.grid(h = -1, v = -1)
  }

  ## if(is.factor(x)) {
##     nint   <- length(levels(x))
##     breaks <- seq(0.5, length = length(levels(x))+1)

##   } else {
##     nint      <- round(log2(length(x))+1)
##     endpoints <- range(x[!is.na(x)])
##     breaks    <- do.breaks(endpoints, nint)
##   }
  # panel.qqmathline(x,breaks,...)

  panel.qqmathline(x,
               lty = abllty,
               lwd = abllwd,
               col.line = ablcol,
                   ...)
  panel.qqmath(x,
               #breaks,
               pch=pch,
               col=col,
               cex=cex,
               lty = abllty,
               lwd = abllwd,
               col.line = ablcol,
               ...)
}

  
