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

addid <- function(x,y,
                  ids   = ids,
                  idsmode=NULL,
                  idsext =0.05,
                  idscex= 0.7,
                  idsdir= "both",
                  gridmode=TRUE
                  )  {

  textfun <- "text"
  if(gridmode) textfun <- "ltext"
  
  if(!is.null(idsmode)) {
    do.call(textfun,list(x,y,ids,cex=idscex))
  } else {
    idd  <- ids
    yres <- residuals(loess(y~x))
    xcut <- cut(x,breaks=4,include.lowest=TRUE)

    ## Determine the number of points to plot
    if(idsext < 1) { # Fraction of total number of points
      idsext    <- ceiling(length(x)*idsext/4)
    }

    for(pp in levels(xcut)) {
      sel   <- xcut == pp
      yyres <- yres[sel]
      xx    <- x[sel]
      yy    <- y[sel]
      iidd  <- idd[sel]
      
      ord     <- order(yyres)
      ordy    <- yy[ord]
      ordx    <- xx[ord]
      ordidd  <- iidd[ord]

      if(!is.null(idsdir)) {
      ## Lower extreme
      if(idsdir=="both" || idsdir=="down") {
        do.call(textfun,list(ordx[1:idsext],
                             ordy[1:idsext],
                             ordidd[1:idsext],
                             cex=idscex))
      }
      ## Upper extreme
      if(idsdir=="both" || idsdir=="up") {
        ll <- length(ordx)
        
        ## Sanity check on ll added by Justin 10/10/05
        ## to prevent subscript errors in conditioned plots
        if (ll!=0) {
        try(
        do.call(textfun,list(ordx[(ll-idsext+1):ll],
                             ordy[(ll-idsext+1):ll],
                             ordidd[(ll-idsext+1):ll],
                             cex=idscex)), silent=T)
        }
        
      }
      }
    }
  }
}

  
