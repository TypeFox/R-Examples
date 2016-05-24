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

"absval.wres.vs.idv" <-
  function(object,
           idv="idv",
           wres="Default",
           ylb  = "Default",
           smooth       = TRUE,
           idsdir       = "up",
           type         = "p",
           ...) {

    if(is.null(check.vars(idv,object,silent=FALSE))) {
      return()
    }

    if(!is.na(match(wres,"Default"))) {
      if(is.null(check.vars("cwres",object,silent=T))) {
        if(is.null(check.vars("wres",object,silent=F))){
          return()
        }else{
          wres="wres"
        }
      } else {
        wres = "cwres"
      }
    } else {
      if(is.null(check.vars(wres,object,silent=FALSE))) {
        return()
      }
    }
    
    
    if(is.null(xvardef(idv,object))){
      xvar <- idv
    }else{
      xvar <- xvardef(idv,object)
    }
    if(is.null(xvardef(wres,object))){
      yvar <- wres
    }else{
      yvar <- xvardef(wres,object)
    }

    if(!is.na(match(ylb,"Default"))) {
      ylb <- paste("|",yvar, "|", sep="")
    } 
    
    xplot <- xpose.plot.default(xvar,
                                yvar,
                                object,
                                ylb=ylb,
                                funy="abs",
                                idsdir=idsdir,
                                smooth=smooth,
                                type = type,
                                ...)
    
    return(xplot)
  }

