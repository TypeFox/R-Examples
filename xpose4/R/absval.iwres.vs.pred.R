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

"absval.iwres.vs.pred" <-
  function(object,
           ylb  = "|IWRES|",
           smooth       = TRUE,
           idsdir       = "up",
           type         = "p",
           ...) {

    if(is.null(xvardef("iwres",object)) ||
       is.null(xvardef("pred",object))) {
      cat("The required variables are not set in the database!\n")
      return()
    }


    xplot <- xpose.plot.default(xvardef("pred",object),
                                xvardef("iwres",object),
                                object,
                                ylb=ylb,
                                funy="abs",
                                idsdir=idsdir,
                                smooth=smooth,
                                type = type,
                                ...)
    
    return(xplot)
  }

