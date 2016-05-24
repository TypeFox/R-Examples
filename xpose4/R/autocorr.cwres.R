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

"autocorr.cwres" <-
  function(object,
           #ylb  = "|WRES|",
           #idsdir="up",
           type="p",
           smooth=TRUE,
           ids=F,
           main = "Default",
           ...) {

    if(is.null(check.vars(c("cwres"),
                          object,silent=FALSE))) {
      return()
    }
    

    default.plot.title <- paste("Autocorrelation of ",xlabel(xvardef("cwres",object),object),
                                sep="")
    plotTitle <- xpose.multiple.plot.title(object=object,
                                           plot.text = default.plot.title,
                                           main=main,
                                           ...)

    xplot <- xpose.plot.default(xvardef("cwres",object),
                                xvardef("cwres",object),
                                object,
                                #ylb=ylb,
                                #idsdir=idsdir,
                                type=type,
                                smooth=smooth,
                                ids=ids,
                                autocorr=TRUE,
                                main=plotTitle,
                                ...)

    return(xplot)
  }
