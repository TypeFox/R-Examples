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

"xpose.panel.splom" <-
  function(x, y, object, subscripts,
           onlyfirst=TRUE,
           inclZeroWRES=FALSE,
           type = "p",
           col  = object@Prefs@Graph.prefs$col,
           pch  = object@Prefs@Graph.prefs$pch,
           cex  = object@Prefs@Graph.prefs$cex,
           lty  = object@Prefs@Graph.prefs$lty,
           lwd  = object@Prefs@Graph.prefs$lwd,
           
           smooth=  TRUE, 
           smlwd = object@Prefs@Graph.prefs$smlwd, 
           smlty = object@Prefs@Graph.prefs$smlty, 
           smcol = object@Prefs@Graph.prefs$smcol, 
           smspan= object@Prefs@Graph.prefs$smspan,
           smdegr= object@Prefs@Graph.prefs$smdegr,
           
           lmline = NULL,
           lmlwd = object@Prefs@Graph.prefs$lmlwd ,
           lmlty = object@Prefs@Graph.prefs$lmlty ,
           lmcol = object@Prefs@Graph.prefs$lmcol ,
           
           grid = object@Prefs@Graph.prefs$grid,
           
           ##scales = list(),
           groups = NULL,
           ... ) {
    
    if(grid != FALSE) {
      panel.grid(h = -1, v = -1)
    }
    
    if(any(is.null(groups))) {
      panel.splom(x, y,
                  col   =col,
                  pch   =pch,
                  lty   =lty,
                  type  =type,
                  cex   = cex,
                  lwd   = lwd,
                  ...
                  )
    } else {
      ord <- order(x)
      panel.superpose(x[ord],
                      y[ord],
                      subscripts[ord],
                      pch   =pch,
                      cex   = cex,
                      lty   =lty,
                      type  =type,
                      lwd   = lwd,
                      groups=groups
                      )
    }
    
    if(!any(is.null(smooth))) {
      if(!is.factor(x) & !is.factor(y)){
        panel.loess(x,y,
                    span  = smspan,
                    degree= smdegr,
                    col   = smcol,
                    lwd   = smlwd,
                    lty   = smlty )
      } else {
        if(is.factor(x) & !is.factor(y)){
          panel.average(x, y,
                        fun = median,
                        horizontal = FALSE,
                        lwd=smlwd, lty=smlty, col=smcol,
                        col.line=smcol,
                        #type="l",
                        ...)
        }
        if(!is.factor(x) & is.factor(y)){
          panel.linejoin(x, y, fun = median, horizontal = TRUE,
                         lwd=smlwd, lty=smlty, col=smcol,
                         col.line=smcol, #type=smlty,
                         ...)
        }
        
      }
       

    }
    if(!any(is.null(lmline))) {
      if(!is.factor(x) & !is.factor(y)){
        panel.abline(lm(y~x),
                     col   = lmcol,
                     lwd   = lmlwd,
                     lty   = lmlty
                     )
      }
    }
  }
		                
