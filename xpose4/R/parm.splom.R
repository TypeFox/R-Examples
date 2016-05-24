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

parm.splom <- function(object,  
                      main = xpose.multiple.plot.title(object=object,
                        plot.text = "Scatterplot matrix of parameters",
                        ...),
                      varnames  = NULL,
                                        #xlb = NULL,
                                        #ylb = NULL,
                      onlyfirst=TRUE,
                                        #inclZeroWRES=FALSE,
                                        #subset=xsubset(object),
                      smooth = TRUE,
                      lmline = NULL,
                                        #groups = NULL,
                                        #main.cex=NULL,
                      ...) {
  
  if(any(is.null(xvardef("parms",object)))) {
    return(cat("Parameters are not defined in the current database!\n"))
  }
  
  if(is.null(varnames)) {
    varnames <- c()
    for (i in xvardef("parms", object)) {
      varnames <- c(varnames, xlabel(i, object))
    }
  }
  
  xplot <- xpose.plot.splom(xvardef("parms", object),
                            object,
                            varnames=varnames,
                            main = main,
                            onlyfirst = onlyfirst,
                                        #inclZeroWRES = inclZeroWRES,
                                        #subset = subset,
                                        #groups = groups,
                            smooth = smooth,
                            lmline = lmline,
                                        #ylb = ylb,
                                        #xlb = xlb,
                            ...)
  
  return(xplot)

}
